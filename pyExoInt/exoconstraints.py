from exomodel import exomodel

# get pandas DataFrame of results
result_df = exomodel('stellar_abu_files/examplestar.txt',
                     Nsim = 1e4,
                     refsolar = 'A21',
                     )

print(result_df)
