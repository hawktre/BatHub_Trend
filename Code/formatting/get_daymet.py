# save as get_daymet.py
import pandas as pd
import pydaymet

# load your site file
sites = pd.read_csv("DataRaw/covariates/daymet/daymet_batch.csv",
                    parse_dates=["night"])

out = []

for _, row in sites.iterrows():
    lat, lon, date = row["latitude"], row["longitude"], row["night"].date()
    df = pydaymet.get_bycoords(
        coords=[(lon, lat)],
        dates=(date, date),
        variables=["tmin","prcp","vp","dayl"]
    )
    # pydaymet returns multi-index: coords and time → extract first
    df = df.reset_index()
    df["site"] = row["site"]
    out.append(df)

result = pd.concat(out, ignore_index=True)
result.to_csv("all_daymet.csv", index=False)
