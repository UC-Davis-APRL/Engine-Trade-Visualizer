import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from Analysis.Evaluation import twoDimensionalSweep

def strip_units(df):

    units = {}
    clean = {}

    for col in df.columns:
        sample = df[col].iloc[0]

        if hasattr(sample, 'units'):  # is a pint quantity
            units[col] = str(sample.units)
            clean[col] = df[col].apply(lambda x: x.magnitude)
        else:
            units[col] = ""
            clean[col] = df[col]

    return pd.DataFrame(clean), units

def df_to_grid(df_clean, param1, param2, color_output):
    pivot = df_clean.pivot(index=param2, columns=param1, values=color_output)
    return pivot


def plot_heatmap(df, param1, param2, color_output, hover_outputs):
    df_clean, units = strip_units(df)

    def label(col):
        u = units.get(col, "")
        return f"{col} ({u})" if u else col

    pivot = df_clean.pivot(index=param2, columns=param1, values=color_output)

    fig = px.imshow(
        pivot.values,
        labels=dict(
            x=label(param1),
            y=label(param2),
            color=label(color_output)
        ),
        x=pivot.index.values,
        y=pivot.columns.values,
        aspect="auto",
        origin="lower",
        title=f"{label(color_output)} vs {label(param1)} & {label(param2)}",
        color_continuous_scale="Viridis"
    )

    # build customdata stack from hover outputs
    hover_arrays = []
    for col in hover_outputs:
        pivot_col = df_clean.pivot(index=param2, columns=param1, values=col).values
        hover_arrays.append(pivot_col)

    if hover_arrays:
        fig.update_traces(
            customdata=np.stack(hover_arrays, axis=-1)
        )

        # build hovertemplate dynamically from whatever outputs the user selected
        hover_template = (
            f"{label(param1)}: %{{x:.3f}}<br>"
            f"{label(param2)}: %{{y:.3f}}<br>"
            f"{label(color_output)}: %{{z:.3f}}<br>"
        )

        formats = {
            "Isp": ".2f",
            "Fuel Volume": ".6e",  # scientific for small values
        }

        for i, col in enumerate(hover_outputs):
            fmt = formats.get(col, ".3f")  # default to .3f if not specified
            hover_template += f"{label(col)}: %{{customdata[{i}]:{fmt}}}<br>"

        hover_template += "<extra></extra>"

        fig.update_traces(hovertemplate=hover_template)

    fig.update_layout(width=800, height=700)

    return fig

st.title("Rocket God Calculator")

with st.sidebar:
    st.header("Initial Configuration")

    FuelList = ['RP-1 (RPL)', 'METHANE', 'PROPANE', 'BUTANE(2,2-BISDIFLUOROAMINO)', 'ISOPROPYL ALCOHOL','HYDROGEN (CRYOGENIC)', 'ETHANOL']
    OxidizerList = ['OXYGEN (LIQUID)', 'OXYGEN (GAS)', 'NITROUS OXIDE', 'AIR (DRY AT SEA LEVEL)', 'CHLORINE', 'HYDROGEN PEROXIDE (100 PC)']
    fuel        = st.selectbox("Choose your Fuel", FuelList, index = 0)
    ox        = st.selectbox("Choose your Oxidizer", OxidizerList, index = 0)

    of        = st.number_input("O/F",          value=2.5)
    pc        = st.number_input("Pc (psia)",    value=300.0)

    thrust    = st.number_input("Thrust (lbf)", value=5000.0)
    burn_time = st.number_input("Burn Time (s)",value=30.0)

    st.header("Sweep Settings")
    params = ["OF", "Pc", "thrust", "burn_time"]
    param1 = st.selectbox("Sweep Param 1", options=params, index=0)
    param2 = st.selectbox("Sweep Param 2", options=params, index=1)

    min1, max1 = st.slider("Range 1", 0.0, 1000.0, (200.0, 400.0))
    min2, max2 = st.slider("Range 2", 0.0, 10.0,   (2.0,   3.0))
    n_points   = st.slider("Points per axis", 5, 50, 20)

    selected_outputs = st.multiselect(
        "Outputs",
        options=["Isp", "Fuel Volume", "Oxidizer Volume"],
        default=["Isp", "Fuel Volume"]
    )

if "df" not in st.session_state:
    st.session_state.df = None

if st.button("Run Sweep"):
    if not selected_outputs:
        st.warning("Select at least one output")
        st.stop()

    config = {
        "OF": 5,
        "Pc": 20,
        "mdot": 2,
        "Thrust": None,
        "burnTime": 10,
        "vehicleRadius": 2,
        "oxName": 'LOX',
        "fuelName":'RP_1'
    }

    range1 = np.linspace(min1, max1, n_points)
    range2 = np.linspace(min2, max2, n_points)

    with st.spinner("Running sweep..."):
        st.session_state.df = twoDimensionalSweep(config, param1, range1, param2, range2, selected_outputs)

if st.session_state.df is not None:
    df = st.session_state.df
    print(df)

    color_output  = st.selectbox("Heatmap Color", options=selected_outputs)
    hover_outputs = [o for o in selected_outputs if o != color_output]

    fig = plot_heatmap(df, param1, param2, color_output, hover_outputs)
    st.plotly_chart(fig, use_container_width=True)
    st.dataframe(df)