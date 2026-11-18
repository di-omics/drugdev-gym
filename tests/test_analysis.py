import numpy as np

from drugdevgym.analysis import analyte_zprime, estimate_half_life, fit_4pl, invert_4pl, zprime
from drugdevgym.chem import Kit, hill, remaining


def test_zprime_matches_the_definition():
    pos, neg = [10.0, 11.0, 12.0], [100.0, 101.0, 102.0]
    assert np.isclose(zprime(pos, neg), 1 - 3 * (1.0 + 1.0) / 90.0)
    assert np.isnan(zprime([1.0, 2.0], [5.0, 6.0, 7.0]))


def test_fit_4pl_recovers_a_rising_curve():
    logx = np.linspace(-1, 4, 10)
    y = 5.0 + 95.0 * hill(10.0 ** logx, 10.0 ** 1.7, 1.2)
    fit = fit_4pl(logx, y)
    assert abs(fit["logec50"] - 1.7) <= 0.02 and abs(fit["hill"] - 1.2) <= 0.1
    assert abs(fit["bottom"] - 5.0) < 1.0 and abs(fit["top"] - 100.0) < 1.0


def test_fit_4pl_handles_a_falling_standard_curve_and_inverts_it():
    kit = Kit(k=30.0)
    std = np.array([1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0])
    fit = fit_4pl(np.log10(std), kit.ratio(std))
    back = invert_4pl(fit, kit.ratio(np.array([5.0, 30.0, 200.0])))
    assert np.allclose(back, [5.0, 30.0, 200.0], rtol=0.1)


def test_fit_4pl_needs_four_points():
    assert fit_4pl([0, 1, 2], [1, 2, 3]) is None


def test_competitive_ratio_falls_as_analyte_rises():
    kit = Kit()
    assert kit.ratio(0.0) == kit.rmax and kit.ratio(1e9) < kit.rmin + 1.0
    assert np.all(np.diff(kit.ratio(np.array([1.0, 10.0, 100.0]))) < 0)


def test_analyte_zprime_sees_saturation_that_raw_zprime_misses():
    kit, rng = Kit(k=30.0), np.random.default_rng(0)
    std = [{"std_nM": s, "ratio": float(kit.ratio(s))} for s in (1, 3, 10, 30, 100, 300, 1000, 3000)]

    def read(level):
        return kit.ratio(level * (1 + rng.normal(0, 0.06, 16))) + rng.normal(0, 30.0, 16)

    good = analyte_zprime(std, read(150.0), read(10.0))
    saturated_pos, saturated_neg = read(3000.0), read(300.0)
    assert good > 0.5
    assert analyte_zprime(std, saturated_pos, saturated_neg) < zprime(saturated_pos, saturated_neg)
    assert analyte_zprime(std, saturated_pos, saturated_neg) < 0.5


def _samples(half_life, times, floor=0.0):
    return [{"id": "c", "replicate": 1, "time_min": float(t), "area": 1000.0 * remaining(half_life, t) + floor}
            for t in times]


def test_half_life_is_recovered_without_noise():
    est = estimate_half_life(_samples(42.0, [0, 15, 30, 60, 120]))
    assert np.isclose(est["c"], 42.0, rtol=0.01)


def test_stable_compound_reports_the_upper_limit():
    assert estimate_half_life(_samples(1e6, [0, 60, 360]))["c"] == 1440.0


def test_lloq_keeps_a_noise_floor_from_flattening_the_slope():
    samples = _samples(8.0, [0, 5, 15, 30, 60, 120, 240, 360], floor=3.0)
    with_lloq = estimate_half_life(samples, use_lloq=True)["c"]
    without = estimate_half_life(samples, use_lloq=False)["c"]
    assert abs(with_lloq - 8.0) < 1.5 and without > 3 * with_lloq


def test_everything_below_lloq_gives_a_bound_not_a_fit():
    est = estimate_half_life(_samples(2.0, [0, 360]))["c"]
    assert est <= 360 / np.log2(50) + 1e-6


def test_flat_data_gives_no_potency():
    from drugdevgym.analysis import estimate_pec50
    rng = np.random.default_rng(0)
    wells = [{"role": "cpd", "id": "c", "conc_nM": 10.0 ** x, "std_nM": None, "ratio": 5000.0 + rng.normal(0, 50)}
             for x in np.linspace(-1, 4, 10)]
    assert estimate_pec50(wells, use_standards=False) == {}
    assert np.isnan(invert_4pl({"bottom": 1.0, "top": 1.0, "logec50": 0.0, "hill": 1.0}, [1.0, 2.0])).all()


def test_half_life_needs_a_time_course_and_ignores_bad_points():
    one_point = [{"id": "c", "replicate": 1, "time_min": 0.0, "area": 1000.0}]
    assert estimate_half_life(one_point) == {}
    bad_t0 = _samples(42.0, [0, 15, 30, 60, 120]) + [{"id": "c", "replicate": 2, "time_min": 0.0, "area": float("nan")}]
    assert np.isclose(estimate_half_life(bad_t0)["c"], 42.0, rtol=0.01)
    zero_is = [{**s, "is_area": 0.0} for s in _samples(42.0, [0, 15, 30])]
    assert estimate_half_life(zero_is) == {}
