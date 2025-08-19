import pandas as pd
df = pd.read_parquet("/home/ilay/hera_unittest_data/measurements/meteorology/highfreqdata/slicedYamim_sonic.parquet")
print(df.columns)
