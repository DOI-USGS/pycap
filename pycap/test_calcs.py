from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import pycap
from pycap.utilities import Q2ts, create_timeseries_template


# homepath = Path(getcwd())
# datapath = homepath / 'tests' / 'data'
datapath = Path("pycap/tests/data")

create_timeseries_template(
    filename=datapath / "test_ts.csv",
    well_ids=[f"well{i}" for i in range(1, 6)],
)


@pytest.fixture
def theis_results():
    excel_file = datapath / "HighCap_Analysis_Worksheet_Example.xlsm"
    p = pd.read_excel(
        excel_file,
        sheet_name="Property_Drawdown_Analysis",
        usecols="C:D",
        skiprows=7,
        index_col=0,
    )
    # read the two Q values and convert to CFD
    Q = [
        float(i) * pycap.GPM2CFD
        for i in [
            p.loc["Pumping Rate Well #1 (gpm)"].iloc[0],
            p.loc["Pumping Rate Well #2 (gpm)"].iloc[0],
        ]
    ]
    S = float(p.loc["Storage Coefficient (unitless)"].iloc[0])
    T = float(p.loc["Transmissivity (ft2/day)"].iloc[0])
    time = float(p.loc["Theis Time of Pumping (days)"].iloc[0])
    params = {"Q": Q, "S": S, "T": T, "time": time}
    theis_res = pd.read_excel(
        excel_file,
        sheet_name="Theis",
        skiprows=9,
        usecols=("A:B,H:I"),
        names=["well1_dd", "well1_r", "well2_dd", "well2_r"],
    )

    return {"params": params, "theis_res": theis_res}


@pytest.fixture
def ward_lough_test_data():
    s1_test = pd.read_csv(datapath / "s1_test.csv", index_col=0)
    s2_test = pd.read_csv(datapath / "s2_test.csv", index_col=0)
    dQ1_test = pd.read_csv(datapath / "dQ1_test.csv", index_col=0)
    dQ2_test = pd.read_csv(datapath / "dQ2_test.csv", index_col=0)
    params = {
        "T1": 100,
        "T2": 100,
        "S1": 1000,
        "S2": 1,
        "width": 1,
        "Q": 125,
        "dist": 100,
        "streambed_thick": 10,
        "streambed_K": 1,
        "aquitard_thick": 10,
        "aquitard_K": 0.01,
        "x": 50,
        "y": 100,
    }
    return {
        "s1_test": s1_test,
        "s2_test": s2_test,
        "dQ1_test": dQ1_test,
        "dQ2_test": dQ2_test,
        "params": params,
    }


@pytest.fixture
def walton_results():
    excel_file = datapath / "HighCap_Analysis_Worksheet_Example.xlsm"
    walton_res = pd.read_excel(
        excel_file,
        sheet_name="Stream#1_Depletion",
        skiprows=104,
        usecols=("C,M:N,R,AB:AC,AK"),
        names=[
            "t_well",
            "dep1",
            "dep2",
            "t_image",
            "rch1",
            "rch2",
            "total_dep",
        ],
    )
    p = pd.read_excel(
        excel_file,
        sheet_name="Stream#1_Depletion",
        skiprows=70,
        nrows=30,
        usecols=("B:D"),
        names=["par", "v1", "v2"],
        index_col=0,
    )
    Q = p.loc["Q - Pumping rate (ft^3/dy)"].T.values
    S = p.loc["S - Storage (unitless)"].T.values
    dist = p.loc["a - Distance (feet)"].T.values
    T_gpd_ft = p.loc["T - Transmissivity (gpd/ft)"].T.values
    # a little trickery to get the index of the time array for start and end of each well
    Q_start_day = [
        pd.to_datetime(i).day_of_year - 1
        for i in p.loc["  First Day of Annual Pumping ="].T.values
    ]
    Q_end_day = [
        pd.to_datetime(i).day_of_year - 1
        for i in p.loc["  Last Day of Annual Pumping ="].T.values
    ]
    params = {
        "Q": Q,
        "S": S,
        "T_gpd_ft": T_gpd_ft,
        "Q_start_day": Q_start_day,
        "Q_end_day": Q_end_day,
        "dist": dist,
    }
    return {"params": params, "walton_res": walton_res}


@pytest.fixture
def project_spreadsheet_results():
    excel_file = datapath / "HighCap_Analysis_Worksheet_Example.xlsm"
    # read in common parameters
    p = pd.read_excel(
        excel_file,
        sheet_name="Property_Drawdown_Analysis",
        skiprows=7,
        nrows=12,
        usecols=("C:D"),
        index_col=0,
    )
    p1 = pd.read_excel(
        excel_file,
        sheet_name="Property_Drawdown_Analysis",
        skiprows=19,
        nrows=12,
        usecols=("C:D"),
        index_col=0,
    )
    p2 = pd.read_excel(
        excel_file,
        sheet_name="Property_Drawdown_Analysis",
        skiprows=31,
        nrows=12,
        usecols=("C:D"),
        index_col=0,
    )
    p3 = pd.read_excel(
        excel_file,
        sheet_name="Property_Drawdown_Analysis",
        skiprows=57,
        nrows=50,
        usecols=("C:F"),
        index_col=0,
    )
    p4 = pd.read_excel(
        excel_file,
        sheet_name="Cumulative_Impact_Analysis",
        skiprows=22,
        nrows=5,
        usecols=("C:D"),
        index_col=0,
    )
    p5 = pd.read_excel(
        excel_file,
        sheet_name="Cumulative_Impact_Analysis",
        skiprows=36,
        nrows=10,
        usecols=("H:I"),
        index_col=0,
    )

    params = {
        "T": p.loc["Transmissivity (ft2/day)"].values[0],
        "S": p.loc["Storage Coefficient (unitless)"].values[0],
        "Q1_gpm": p.loc["Pumping Rate Well #1 (gpm)"].values[0],
        "Q2_gpm": p.loc["Pumping Rate Well #2 (gpm)"].values[0],
        "w1muni_dist": p3.loc["Distance from Well #1 to Municpal Well"].values[
            0
        ],
        "w2muni_dist": p3.loc["Distance from Well #2 to Municpal Well"].values[
            0
        ],
        "w1sprng1_dist": p3.loc["Distance from Well #1 to Spring"].values[0],
        "w2sprng1_dist": p3.loc["Distance from Well #2 to Spring"].values[0],
        "muni_dd_combined_proposed": p3.loc[
            "Distance from Well #1 to Municpal Well"
        ].values[-1],
        "sprng1_dd_combined_proposed": p3.loc[
            "Distance from Well #1 to Spring"
        ].values[-1],
        "well1_5ftdd_loc": p3.loc[" Well #1 5-ft Drawdown (feet)"].values[0],
        "well1_1ftdd_loc": p3.loc[" Well #1 1-ft Drawdown (feet)"].values[0],
        "theis_p_time": p.loc["Theis Time of Pumping (days)"].values[0],
        "stream_name_1": p1.loc["Stream Name"].values[0],
        "stream_name_2": p2.loc["Stream Name"].values[0],
        "depl_pump_time": p1.loc[
            "Stream Depletion Duration Period (Days)"
        ].values[0],
        "w1s1_dist": p1.loc["Well #1 - Distance to Stream (feet)"].values[0],
        "w1s1_appor": p1.loc[
            "Well #1 - Fraction Intercepting Stream (.1-1)"
        ].values[0],
        "w2s1_dist": p1.loc["Well #2 - Distance to Stream (feet)"].values[0],
        "w2s1_appor": p1.loc[
            "Well #2 - Fraction Intercepting Stream (.1-1)"
        ].values[0],
        "w1s2_dist": p2.loc["Well #1 - Distance to Stream (feet)"].values[0],
        "w1s2_appor": p2.loc[
            "Well #1 - Fraction Intercepting Stream (.1-1)"
        ].values[0],
        "w2s2_dist": p2.loc["Well #2 - Distance to Stream (feet)"].values[0],
        "w2s2_appor": p2.loc[
            "Well #2 - Fraction Intercepting Stream (.1-1)"
        ].values[0],
        "s1_4yr_depl_cfs": p3.loc[
            "Stream #1 depletion after year 4 (cfs)"
        ].values[0],
        "s2_4yr_depl_cfs": p3.loc[
            "Stream #2 depletion after year 4  (cfs)"
        ].values[0],
        "muni_dd_total_combined": p4.loc[
            "Cumulative Impact Drawdown (ft)"
        ].values[0],
        "stream1_depl_existing": p5.iloc[0].values[0],
        "stream1_depl_total_combined": p5.iloc[3].values[0],
    }
    return params


@pytest.fixture
def hunt_03_results():
    # read in results from STRMDEPL08 example run
    flname = datapath / "example03.plt"
    strmdepl08_df = pd.read_csv(flname, sep=r"\s+")
    strmdepl08_df.index = (
        strmdepl08_df.index + 1
    )  # adjust index to match python output
    strmdepl08_df["ratio08"] = strmdepl08_df["QS"] / strmdepl08_df["QWELL"]
    time = [50, 100, 200, 300]
    checkvals = [strmdepl08_df.loc[x]["ratio08"] for x in time]
    return {"time": time, "checkvals": checkvals}


def test_project_spreadsheet(project_spreadsheet_results):
    import pycap
    from pycap.wells import Well

    pars = project_spreadsheet_results
    # set up the Project with multiple wells and multiple streams and make calculations
    well1 = Well(
        T=pars["T"],
        S=pars["S"],
        Q=Q2ts(pars["depl_pump_time"], 5, pars["Q1_gpm"]) * pycap.GPM2CFD,
        depletion_years=5,
        theis_dd_days=pars["theis_p_time"],
        depl_pump_time=pars["depl_pump_time"],
        stream_dist={
            pars["stream_name_1"]: pars["w1s1_dist"],
            pars["stream_name_2"]: pars["w1s2_dist"],
        },
        drawdown_dist={"muni": pars["w1muni_dist"]},
        stream_apportionment={
            pars["stream_name_1"]: pars["w1s1_appor"],
            pars["stream_name_2"]: pars["w1s2_appor"],
        },
    )
    well2 = Well(
        T=pars["T"],
        S=pars["S"],
        Q=Q2ts(pars["depl_pump_time"], 5, pars["Q2_gpm"]) * pycap.GPM2CFD,
        depletion_years=5,
        theis_dd_days=pars["theis_p_time"],
        depl_pump_time=pars["depl_pump_time"],
        stream_dist={
            pars["stream_name_1"]: pars["w2s1_dist"],
            pars["stream_name_2"]: pars["w2s2_dist"],
        },
        drawdown_dist={"muni": pars["w2muni_dist"]},
        stream_apportionment={
            pars["stream_name_1"]: pars["w2s1_appor"],
            pars["stream_name_2"]: pars["w2s2_appor"],
        },
    )
    dd1 = well1.drawdown["muni"][well1.theis_dd_days]
    dd2 = well2.drawdown["muni"][well2.theis_dd_days]

    assert np.allclose(dd1 + dd2, pars["muni_dd_combined_proposed"], atol=0.1)

    depl1 = well1.depletion
    depl1 = {k: v / 3600 / 24 for k, v in depl1.items()}
    depl2 = well2.depletion
    depl2 = {k: v / 3600 / 24 for k, v in depl2.items()}
    stream1_max_depl = np.max(depl1[pars["stream_name_1"]]) + np.max(
        depl2[pars["stream_name_1"]]
    )
    stream2_max_depl = np.max(depl1[pars["stream_name_2"]]) + np.max(
        depl2[pars["stream_name_2"]]
    )
    assert np.allclose(stream1_max_depl, pars["s1_4yr_depl_cfs"], atol=1e-2)
    assert np.allclose(stream2_max_depl, pars["s2_4yr_depl_cfs"], atol=1e-2)


def test_theis(theis_results):
    """Test for the theis calculations - compared with two wells at multiple distances
        in the example spreadsheet

    Args:
        theis_results (@fixture, dict): parameters and results from example spreadsheet
    """

    pars = theis_results["params"]
    dist = theis_results["theis_res"].well1_r

    time = pars["time"]
    dd = [
        pycap.theis_drawdown(pars["T"], pars["S"], time, dist, currQ)
        for currQ in pars["Q"]
    ]
    assert np.allclose(dd[0], theis_results["theis_res"].well1_dd, atol=0.5)
    assert np.allclose(dd[1], theis_results["theis_res"].well2_dd, atol=0.7)


def test_distance():
    from pycap import analysis_project as ap

    assert np.isclose(
        ap._loc_to_dist([89.38323, 43.07476], [89.38492, 43.07479]),
        450.09,
        atol=0.1,
    )
    #  ([2,3],[9,32.9]), 30.70846788753877)


def test_glover_depletion():
    """Test for the glover calculations
    against the Glover & Balmer (1954) paper
    """
    dist = [1000, 5000, 10000]
    Q = 1
    time = 365 * 5  # paper evaluates at 5 years in days
    K = 0.001  # ft/sec
    D = 100  # thickness in feet
    T = K * D * 24 * 60 * 60  # converting to ft/day
    S = 0.2
    Qs = pycap.glover_depletion(T, S, time, dist, Q)
    assert not any(np.isnan(Qs))
    assert np.allclose(Qs, [0.9365, 0.6906, 0.4259], atol=1e-3)


def test_sdf():
    """Test for streamflow depletion factor
    using values from original Jenkins (1968) paper
    https://doi.org/10.1111/j.1745-6584.1968.tb01641.x
    note Jenkins rounded to nearest 10 (page 42)
    """

    dist = 5280.0 / 2.0
    T = 5.0e4 / 7.48
    S = 0.5
    sdf = pycap.sdf(T, S, dist)
    assert np.allclose(sdf, 520, atol=1.5)


# def test_well():
#     from pycap import wells
#     w = wells.Well('pending',
#                    )


def test_walton_depletion(walton_results):
    """Test of a single year to be sure the Walton calculations are made correctly

    Args:
        walton_results ([type]): [description]
    """
    res = walton_results["walton_res"]
    pars = walton_results["params"]
    # prepend a zero time to be sure that logic is tested
    cols = res.columns
    res = pd.DataFrame(
        np.insert(res.values, 0, values=[0] * len(res.columns), axis=0)
    )
    res.columns = cols
    dep = {}
    rch = {}
    for idx in [0, 1]:
        dep[idx] = pycap.walton_depletion(
            pars["T_gpd_ft"][idx],
            pars["S"][idx],
            [0] + res.t_well,
            pars["dist"][idx],
            pars["Q"][idx],
        )
        rch[idx] = pycap.walton_depletion(
            pars["T_gpd_ft"][idx],
            pars["S"][idx],
            [0] + res.t_image,
            pars["dist"][idx],
            pars["Q"][idx],
        )
    dep_tot = dep[0] - rch[0] + dep[1] - rch[1]
    assert np.allclose(dep[0] / 3600 / 24, res.dep1)
    assert np.allclose(dep[1] / 3600 / 24, res.dep2)
    assert np.allclose(rch[0] / 3600 / 24, -res.rch1)
    assert np.allclose(rch[1] / 3600 / 24, -res.rch2)
    assert np.allclose(dep_tot / 3600 / 24, res.total_dep)


def test_yaml_parsing(project_spreadsheet_results):
    pars = project_spreadsheet_results
    from pycap.analysis_project import Project

    ap = Project(datapath / "example.yml")
    # ap.populate_from_yaml(datapath / 'example.yml')
    # verify that the created well objects are populated with the same values as in the YML file
    assert (
        set(ap.wells.keys()).difference(
            set(["new1", "new2", "Existing_CAFO", "Existing_Irrig"])
        )
        == set()
    )
    assert (
        set(ap._Project__stream_responses.keys()).difference(
            set(["Upp Creek", "no paddle"])
        )
        == set()
    )
    assert (
        set(ap._Project__dd_responses.keys()).difference(
            set(["Muni1", "Sprng1"])
        )
        == set()
    )

    # spot check some numbers
    assert ap.wells["new1"].T == 35
    assert np.isclose(pycap.GPM2CFD * 1000, ap.wells["new2"].Q.iloc[0])
    assert ap.wells["new2"].stream_apportionment["Upp Creek"] == 0.3

    ap.report_responses()

    ap.write_responses_csv()

    agg_results = pd.read_csv(ap.csv_output_filename, index_col=0)
    # read in the CSV file and spot check against the spreadsheet output
    assert np.isclose(
        pars["muni_dd_combined_proposed"],
        agg_results.loc["total_proposed", "Muni1:dd (ft)"],
        atol=0.1,
    )
    assert np.isclose(
        pars["sprng1_dd_combined_proposed"],
        agg_results.loc["total_proposed", "Sprng1:dd (ft)"],
        atol=0.002,
    )
    assert np.isclose(
        pars["stream1_depl_existing"],
        agg_results.loc["total_existing", "Upp Creek:depl (cfs)"],
        atol=0.005,
    )
    assert np.isclose(
        pars["stream1_depl_total_combined"],
        agg_results.loc["total_combined", "Upp Creek:depl (cfs)"],
        atol=0.01,
    )


def test_complex_yml():
    from pycap.analysis_project import Project

    ap = Project(datapath / "example2.yml")
    ap.report_responses()
    ap.write_responses_csv()

    df_ts = pd.read_csv(ap.csv_stream_output_ts_filename, index_col=0)
    df_agg = pd.read_csv(ap.csv_stream_output_filename, index_col=0)

    df_ts_max = df_ts.max().to_frame()
    df_ts_max.rename(columns={0: "raw"}, inplace=True)
    s_cols_exist = [
        i for i in df_ts.columns if ("Spring" in i) & ("93444" not in i)
    ]
    s_cols_prop = [
        i for i in df_ts.columns if ("Spring" in i) & ("93444" in i)
    ]

    e_cols_exist = [
        i for i in df_ts.columns if ("EBranch" in i) & ("93444" not in i)
    ]
    e_cols_prop = [
        i for i in df_ts.columns if ("EBranch" in i) & ("93444" in i)
    ]

    s_cols_tot = s_cols_exist + s_cols_prop
    e_cols_tot = e_cols_exist + e_cols_prop

    df_ts_max["read"] = [
        df_agg.loc[i.split(":")[1], i.split(":")[0]] for i in df_ts_max.index
    ]
    assert all(np.isclose(df_ts_max.raw, df_ts_max["read"]))

    keys = (
        "SpringBrook:proposed",
        "SpringBrook:existing",
        "SpringBrook:combined",
        "EBranchEauClaire:proposed",
        "EBranchEauClaire:existing",
        "EBranchEauClaire:combined",
    )
    vals = (
        s_cols_prop,
        s_cols_exist,
        s_cols_tot,
        e_cols_prop,
        e_cols_exist,
        e_cols_tot,
    )
    for k, v in zip(keys, vals):
        df_agg_val = df_agg.loc[f"total_{k.split(':')[1]}", k.split(":")[0]]
        calc_val = np.max(df_ts[v].sum(axis=1))
        assert np.isclose(df_agg_val, calc_val)

    print("stoked")


def test_run_yml_example():
    from pycap.analysis_project import Project

    yml_file = "example.yml"
    ap = Project(datapath / yml_file)
    ap.report_responses()
    ap.write_responses_csv()


def test_hunt_99_depletion_results_multiple_times():
    """Test of hunt_99_depletion() function in the
    well.py module.  Compares computed stream depletion
    to results from Jenkins (1968) Table 1 and the
    strmdepl08 appendix for dist=1000 and multiple times
    """
    dist = [1000]
    Q = 1
    time = [0, 365 * 5]  # paper evaluates at 5 years in days
    K = 0.001  # ft/sec
    D = 100  # thickness in feet
    T = K * D * pycap.SEC2DAY  # converting to ft/day
    S = 0.2
    rlambda = (
        10000.0  # large lambda value should return Glover and Balmer solution
    )
    # see test_glover for these values.
    Qs = pycap.hunt_99_depletion(
        T, S, time, dist, Q, streambed_conductance=rlambda
    )
    assert not np.atleast_1d(np.isnan(Qs)).any()
    assert np.allclose(Qs, [0, 0.9365], atol=1e-3)

    # check some values with varying time, using t/sdf, q/Q table
    # from Jenkins (1968) - Table 1
    dist = 1000.0
    sdf = dist**2 * S / T
    time = [sdf * 1.0, sdf * 2.0, sdf * 6.0]
    obs = [0.480, 0.617, 0.773]
    Qs = pycap.hunt_99_depletion(
        T, S, time, dist, Q, streambed_conductance=rlambda
    )
    assert not any(np.isnan(Qs))
    assert np.allclose(Qs, obs, atol=5e-3)


def test_hunt_99_depletion_results():
    """Test of hunt_99_depletion() function in the
    well.py module.  Compares computedstream depletion
    to results from Jenkins (1968) Table 1 and the
    strmdepl08 appendix.
    """
    dist = [1000, 5000, 10000]
    Q = 1
    time = 365 * 5  # paper evaluates at 5 years in days
    K = 0.001  # ft/sec
    D = 100  # thickness in feet
    T = K * D * pycap.SEC2DAY  # converting to ft/day
    S = 0.2
    rlambda = (
        10000.0  # large lambda value should return Glover and Balmer solution
    )
    # see test_glover for these values.
    Qs = pycap.hunt_99_depletion(
        T, S, time, dist, Q, streambed_conductance=rlambda
    )
    assert not np.atleast_1d(np.isnan(Qs)).any()
    assert np.allclose(Qs, [0.9365, 0.6906, 0.4259], atol=1e-3)

    # check some values with varying time, using t/sdf, q/Q table
    # from Jenkins (1968) - Table 1
    dist = 1000.0
    sdf = dist**2 * S / T
    time = [sdf * 1.0, sdf * 2.0, sdf * 6.0]
    obs = [0.480, 0.617, 0.773]
    Qs = pycap.hunt_99_depletion(
        T, S, time, dist, Q, streambed_conductance=rlambda
    )
    assert not any(np.isnan(Qs))
    assert np.allclose(Qs, obs, atol=5e-3)

    # Check with lower streambed conductance using
    # values from 28 days of pumping from STRMDEPL08 appendix
    # T = 1,000 ft2/d, L = 100 ft, S = 20 ft/d, d = 500 ft, S = 0.1, and Qw = 0.557 ft3/s (250 gal/min).

    dist = 500.0  # feet
    T = 0.116e-1 * pycap.SEC2DAY  # ft^2/sec to ft^2/day
    S = 0.1
    Q = 0.557 * pycap.SEC2DAY  # cfs to cfd
    time = [10.0, 20.0, 28.0]  # days
    rlambda = 0.231e-03 * pycap.SEC2DAY  # ft/sec to ft/day
    obs = np.array([0.1055, 0.1942, 0.2378]) / 0.5570
    Qs = (
        pycap.hunt_99_depletion(
            T, S, time, dist, Q, streambed_conductance=rlambda
        )
        / Q
    )  # normalize results
    assert not any(np.isnan(Qs))
    assert np.allclose(Qs, obs, atol=5e-3)


def test_hunt_03_depletion_results(hunt_03_results):
    """Test of hunt_03_depletion() function in the
    well.py module.  Compares computed stream depletion
    to results from STRMDEPL08 Fortran code using
    example03.dat as input and producing example03.plt
    """

    dist = 500.0
    T = 0.0115740740740741 * 60.0 * 60.0 * 24.0
    S = 0.001
    Qw = 0.557 * 60 * 60 * 24
    Bprime = 20
    Bdouble = 15
    Kprime = 1.1574074074074073e-05 * 60.0 * 60.0 * 24.0
    sigma = 0.1
    width = 5

    time = np.array([0.0] + hunt_03_results["time"])
    rlambda = Kprime * (width / Bdouble)

    Qs = pycap.hunt_03_depletion(
        T,
        S,
        time,
        dist,
        Qw,
        Bprime=Bprime,
        Bdouble=Bdouble,
        aquitard_K=Kprime,
        sigma=sigma,
        width=width,
        streambed_conductance=rlambda,
    )
    ratios = Qs / Qw

    tol = 0.002  # relative tolerance = 0.2 percent
    res = np.array([0.0] + hunt_03_results["checkvals"])
    res
    np.testing.assert_allclose(ratios, res, rtol=tol)


@pytest.mark.xfail
def test_yml_ts_parsing1():
    from pycap.analysis_project import Project

    # this should fail on the integrity tests
    Project(datapath / "example3.yml")


@pytest.fixture
def SIR2009_5003_Table2_Batch_results():
    """The batch column from Table 2, SIR 2009-5003,
    with the groundwater component of the MI water
    withdrawal screening tool.  This table has
    catchments, distances, apportionment (percent),
    analytical solution, and percent*analytical
    solution.  The analytical solution is computed
    using Hunt (1999)'

    """
    check_df = pd.read_csv(
        datapath / "SIR2009_5003_Table2_Batch.csv", dtype=float
    )
    check_df.set_index("Valley_segment", inplace=True)

    return check_df


def test_WellClass(SIR2009_5003_Table2_Batch_results):
    """Test the Well Class ability to distribute
    depletion using the Hunt (1999) solution and
    inverse distance weighting by comparing the results
    to Table 2 from the SIR 2009-5003.  For the test, the
    distances to the streams and well characteristics
    are provided and passed to the Well object.
    Drawdown and depletion are attributes of the object.

    """
    check_df = SIR2009_5003_Table2_Batch_results
    stream_table = pd.DataFrame(
        (
            {"id": 8, "distance": 14802},
            {"id": 9, "distance": 12609.2},
            {"id": 11, "distance": 15750.5},
            {"id": 27, "distance": 22567.6},
            {"id": 9741, "distance": 27565.2},
            {"id": 10532, "distance": 33059.5},
            {"id": 11967, "distance": 14846.3},
            {"id": 12515, "distance": 17042.55},
            {"id": 12573, "distance": 11959.5},
            {"id": 12941, "distance": 19070.8},
            {"id": 13925, "distance": 10028.9},
        )
    )

    # use inverse-distnace weighting apportionment
    invers = np.array([1 / x for x in stream_table["distance"]])
    stream_table["apportionment"] = (1.0 / stream_table["distance"]) / np.sum(
        invers
    )

    # other properties, for the SIR example streambed conductance was
    # taken from the catchment containing the well
    T = 7211.0  # ft^2/day
    S = 0.01
    Q = 70  # 70 gpm
    stream_table["conductance"] = 7.11855
    pumpdays = int(5.0 * 365)

    # Well class needs a Pandas series for pumping, and units should be
    # cubic feet per day
    Q = pycap.Q2ts(pumpdays, 5, Q) * pycap.GPM2CFD

    # Well class needs dictionaries of properties keyed by the well names/ids
    distances = dict(zip(stream_table.id.values, stream_table.distance.values))
    apportion = dict(
        zip(stream_table.id.values, stream_table.apportionment.values)
    )
    cond = dict(zip(stream_table.id.values, stream_table.conductance.values))

    # make a Well object, specify depletion method
    test_well = pycap.Well(
        T=T,
        S=S,
        Q=Q,
        depletion_years=5,
        depl_method="hunt_99_depletion",
        streambed_conductance=cond,
        stream_dist=distances,
        stream_apportionment=apportion,
    )

    # get depletion
    stream_depl = pd.DataFrame(test_well.depletion)

    # convert to GPM to compare with Table 2 and check
    stream_depl = stream_depl * pycap.CFD2GPM

    five_year = pd.DataFrame(stream_depl.loc[1824].T)
    five_year.rename(columns={1824: "Depletion"}, inplace=True)

    tol = 0.01
    np.testing.assert_allclose(
        five_year["Depletion"].values,
        check_df["Estimated_removal_gpm"].values,
        atol=tol,
    )


def test_hunt_continuous():
    # read in the pumping timeseries and the depletion results included as a column
    flname = datapath / "hunt_test_ts.csv"
    assert flname.exists()
    df = pd.read_csv(flname, index_col=3)
    from pycap.analysis_project import Project

    # only one well in the
    ap = Project(datapath / "hunt_example.yml")

    ap.report_responses()

    ap.write_responses_csv()

    agg_results = pd.read_csv(ap.csv_output_filename, index_col=0)
    # read in the CSV file and check against STRMDEPL08 Appendix 1 output (OFR2008-1166)
    assert np.isclose(
        df.resp_testing.max(),
        agg_results.loc["well1: proposed", "testriver:depl (cfs)"],
        atol=0.001,
    )
    assert np.allclose(
        df.resp_testing.values,
        ap.wells["well1"].depletion["testriver"] / 3600 / 24,
        atol=0.001,
    )


def test_hunt_99_drawdown():
    """Test of hunt_99_drawdown() function in the
    well.py module.
    """
    Q = 1
    dist = 200.0
    T = 1000.0
    S = 0.1
    time = [0.0, 28.0]

    # test if stream conductance is zero
    rlambda = 0
    x = 50.0
    y = 0.0

    ddwn = pycap.hunt_99_drawdown(
        T, S, time, dist, Q, streambed_conductance=rlambda, x=x, y=y
    )
    no_stream = pycap.theis_drawdown(T, S, time, (dist - x), Q)
    assert np.allclose(ddwn, no_stream)


def test_transient_dd():
    # read in the pumping timeseries and the depletion results included as a column
    flname = datapath / "transient_dd_ts.csv"
    assert flname.exists()
    from pycap.analysis_project import Project

    # only one well in the
    ap = Project(datapath / "transient_drawdown.yml")

    ap.report_responses()

    ap.write_responses_csv()

    pd.read_csv(ap.csv_output_filename, index_col=0)


def test_ward_lough_depletion(ward_lough_test_data):
    # note: the parameters defined below are intended to result in the nondimensional
    # parameters corresponding with Fig. 6 in DOI: 10.1061/ (ASCE)HE.1943-5584.0000382.
    allpars = ward_lough_test_data["params"]
    allpars["aquitard_thick"] = 1
    dQ1_test = ward_lough_test_data["dQ1_test"]
    dQ2_test = ward_lough_test_data["dQ2_test"]
    allpars["t"] = dQ2_test.index * 100
    dQ2_test["mod"] = pycap.ward_lough_depletion(**allpars)
    allpars["t"] = dQ1_test.index * 100
    allpars["T1"] = 0.01
    allpars.pop("x")
    allpars.pop("y")
    allpars["aquitard_K"] = 0.001
    dQ1_test["mod"] = pycap.ward_lough_depletion(**allpars)
    assert np.allclose(
        dQ1_test["mod"] / allpars["Q"], dQ1_test["dQ"], atol=0.1
    )

    assert np.allclose(
        dQ2_test["mod"] / allpars["Q"], dQ2_test["dQ"], atol=0.1
    )


def test_ward_lough_drawdown(ward_lough_test_data):
    # note: the parameters defined below are intended to result in the nondimensional
    # parameters corresponding with Fig. 3 in DOI: 10.1061/ (ASCE)HE.1943-5584.0000382.
    allpars = ward_lough_test_data["params"]
    s1_test = ward_lough_test_data["s1_test"]
    s2_test = ward_lough_test_data["s2_test"]
    allpars["t"] = s1_test.index * 100
    s1_test["mod"] = pycap.ward_lough_drawdown(**allpars)[:, 0]
    allpars["t"] = s2_test.index * 100
    s2_test["mod"] = pycap.ward_lough_drawdown(**allpars)[:, 1]
    assert np.allclose(
        s1_test["mod"] * allpars["T2"] / allpars["Q"], s1_test["s"], atol=0.035
    )
    assert np.allclose(
        s2_test["mod"] * allpars["T2"] / allpars["Q"], s2_test["s"], atol=0.035
    )


@pytest.mark.xfail
def test_custom_exception():
    from pycap import theis_drawdown

    # this should raise an exception
    theis_drawdown(1, 1, [1, 2], [1, 2], 5)


def test_complex_well(ward_lough_test_data):
    import pycap
    from pycap import ward_lough_depletion
    from pycap.wells import Well

    # get the test parameters
    allpars = ward_lough_test_data["params"]
    # now run the base solutions for comparisons
    allpars["t"] = list(range(365))
    dep1 = ward_lough_depletion(**allpars)

    # now configure for running through Well object
    allpars["T"] = allpars["T1"]
    allpars["S"] = allpars["S1"]
    allpars["stream_dist"] = None
    allpars["drawdown_dist"] = {"dd1": allpars["dist"]}
    allpars["stream_dist"] = {"resp1": allpars["dist"]}
    allpars["stream_apportionment"] = {"resp1": 1.0}
    allpars["Q"] = pycap.Q2ts(365, 1, allpars["Q"])
    allpars.pop("T1")
    allpars.pop("S1")
    allpars.pop("t")
    allpars.pop("dist")

    w = Well(
        "newwell",
        depl_method="ward_lough_depletion",
        **allpars,
    )
    # athens test - just making sure it runs
    depl = w.depletion
    assert len(depl) > 0

    maxdep = w.max_depletion
    assert len(maxdep) == 1

    # now check against non-Well-object calcs only valid for depletion
    assert np.allclose(dep1[1:], depl["resp1"][1:])
