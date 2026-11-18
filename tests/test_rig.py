from dataclasses import replace

import numpy as np

from drugdevgym.rig import EDGE, Plate, Rig, RigProfile, usable_wells, well_name

EXACT = replace(RigProfile(), cv_floor=0.0, cv_k=0.0, bias_k=0.0, carryover=0.0, dead_volume_ul=0.0)


def rig(profile=EXACT, seed=0, **kw):
    return Rig(profile, np.random.default_rng(seed), **kw)


def filled(volume=100.0, conc=10.0):
    p = Plate("src")
    p.vol[0, 0] = volume
    p.a("x")[0, 0] = conc * volume
    return p


def test_layout_helpers():
    assert well_name((0, 0)) == "A1" and well_name((7, 11)) == "H12"
    assert len(usable_wells(False)) == 96 and len(usable_wells(True)) == 60
    assert not any(EDGE[w] for w in usable_wells(True))


def test_transfer_conserves_volume_and_mass():
    src, dst = filled(), Plate("dst")
    rig().transfer(src, (0, 0), dst, (0, 0), 25.0)
    assert np.isclose(src.vol[0, 0] + dst.vol[0, 0], 100.0)
    assert np.isclose(src.a("x")[0, 0] + dst.a("x")[0, 0], 1000.0)
    assert np.isclose(dst.conc("x", (0, 0)), 10.0)


def test_dead_volume_gives_a_silent_short_aspirate():
    src, dst = filled(volume=20.0), Plate("dst")
    r = rig(replace(EXACT, dead_volume_ul=15.0))
    r.transfer(src, (0, 0), dst, (0, 0), 10.0)
    assert np.isclose(dst.vol[0, 0], 5.0)
    assert r.short_aspirates == 1


def test_unmixed_well_under_delivers_solute():
    got = []
    for seed in range(50):
        src, dst = filled(), Plate("dst")
        src.mixed[0, 0] = False
        rig(seed=seed).transfer(src, (0, 0), dst, (0, 0), 10.0)
        got.append(dst.a("x")[0, 0])
        assert np.isclose(src.a("x")[0, 0] + dst.a("x")[0, 0], 1000.0)  # the rest stays behind
    assert 60.0 < np.mean(got) < 80.0  # 100 if mixed; recovery is about 0.7


def test_carryover_needs_a_reused_tip_and_scales_with_stickiness():
    def second_well(new_tip, sticky):
        src, clean, dst = filled(), Plate("clean"), Plate("dst")
        clean.vol[0, 0] = 100.0
        r = rig(replace(EXACT, carryover=0.01), sticky={"x": sticky})
        r.transfer(src, (0, 0), dst, (0, 0), 10.0)
        r.transfer(clean, (0, 0), dst, (0, 1), 10.0, new_tip=new_tip)
        return dst.a("x")[0, 1]

    assert second_well(True, 1.0) == 0.0
    assert np.isclose(second_well(False, 1.0), 1.0)
    assert np.isclose(second_well(False, 5.0), 5.0)


def test_evaporation_is_worse_on_edges_and_a_lid_helps():
    def after(lid):
        p = Plate("p")
        p.vol[:] = 50.0
        rig(RigProfile()).incubate(p, 120.0, lid=lid)
        return p.vol

    open_plate, lidded = after(False), after(True)
    assert open_plate[EDGE].mean() < open_plate[~EDGE].mean() < 50.0
    assert lidded[EDGE].mean() > open_plate[EDGE].mean()


def test_incubation_mixes_by_diffusion():
    p = Plate("p")
    p.vol[:] = 50.0
    p.mixed[:] = False
    r = rig(RigProfile())
    r.incubate(p, 5.0)
    assert not p.mixed.any()
    r.incubate(p, 20.0)
    assert p.mixed.all()


def test_from_qc_recovers_the_pipetting_terms():
    truth = replace(RigProfile(), cv_floor=0.015, cv_k=0.06, bias_k=-0.05)
    v = np.array([1.0, 2.0, 5.0, 10.0, 50.0, 200.0])
    fit = RigProfile.from_qc(v, [truth.pipette_cv(x) for x in v], [truth.pipette_bias(x) for x in v])
    assert np.isclose(fit.cv_floor, 0.015) and np.isclose(fit.cv_k, 0.06) and np.isclose(fit.bias_k, -0.05)
    assert fit.name == "calibrated" and fit.dead_volume_ul == RigProfile().dead_volume_ul


def test_randomized_rigs_differ_and_are_seeded():
    a = RigProfile.randomized(np.random.default_rng(1))
    b = RigProfile.randomized(np.random.default_rng(1))
    c = RigProfile.randomized(np.random.default_rng(2))
    assert a == b and a != c and a.name == "randomized"
