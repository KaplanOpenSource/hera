# measurements/meteorology/highfreqdata/__main__.py

import pandas as pd

def load_and_describe_data():
    sonic_path = "/home/ilay/hera_unittest_data/measurements/meteorology/highfreqdata/slicedYamim_sonic.parquet"
    trh_path = "/home/ilay/hera_unittest_data/measurements/meteorology/highfreqdata/slicedYamim_TRH.parquet"

    try:
        df_sonic = pd.read_parquet(sonic_path)
        df_trh = pd.read_parquet(trh_path)

        # Safely convert time columns if exist
        for df in [df_sonic, df_trh]:
            for col in ["Time", "timestamp"]:
                if col in df.columns:
                    df[col] = pd.to_datetime(df[col], errors="coerce")  # <- This avoids crashes

        print("🌀 Sonic data sample:")
        print(df_sonic.head(), end="\n\n")

        print("🌡️ TRH data sample:")
        print(df_trh.head(), end="\n\n")

        print("📊 Sonic summary:")
        print(df_sonic.describe(include="all"), end="\n\n")

        print("📊 TRH summary:")
        print(df_trh.describe(include="all"))

    except Exception as e:
        print(f"❌ Error loading or processing data: {e}")

if __name__ == "__main__":
    load_and_describe_data()
