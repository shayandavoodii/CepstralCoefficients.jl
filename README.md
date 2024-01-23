# TimeSeries Cepstral Clustering

This repo contains the Julia implementation of a study, entitled "Distance Measures for Effective Clustering of ARIMA Time-Series"[[1](https://doi.org/10.1109/ICDM.2001.989529)]. The research proposes a approach regarding classifying several time series regarding their inner patterns. Authors used Cepstral Analysis in order to find the patterns underneath the time series.

## How to use?

To use the implementation, you should follow the subsequent steps:

1. Clone the repository using:  

    ```raw
    git clone https://github.com/shayandavoodii/TimeSeries-Cepstral-Clustering.git
    ```

2. Open a Julia session in the "TimeSeries-Cepstral-Clustering" directory and run the following commands:

    ```julia
    using Pkg; Pkg.activate(".")
    include("src/CepstralClustering.jl")
    using .CepstralClustering
    ```

Afterward, you should use the [`cc`](https://github.com/shayandavoodii/TimeSeries-Cepstral-Clustering/blob/main/src/CepstralClustering.jl#L12-L46) function as the main function of the implementation. The aforementioned function calculates `n` number of cepstral coefficients according to the given time series.

## Example

In terms of testing the implementation, a sequence of time series is required. In this case, I use the following assets: `["MSFT", "AAPL", "GOOG", "AMZN", "FB", "TSLA", "NVDA", "INTC", "CSCO", "ADBE"]`. I calculate the first 5 cepstral coefficients to be calculated and compared. I use `p=3` as the lag of $AR$ process.

### Fetch data

```julia
using YFinance
tickers = ["MSFT", "AAPL", "GOOG", "AMZN", "META", "TSLA", "NVDA", "INTC", "CSCO", "ADBE"]
querry = [get_prices(ticker, startdt="2019-01-01", enddt="2020-01-01")["adjclose"] for ticker in tickers]
prices = stack(querry, dims=1);
```

Afterward, I use the [`cc`](https://github.com/shayandavoodii/TimeSeries-Cepstral-Clustering/blob/main/src/CepstralClustering.jl#L12-L46) function to calculate the cepstral coefficients.

### Calculate cepstral coefficients

```julia
n=5
p=3
cepscoefs = cc(prices, p, n)
# 5×10 Matrix{Float64}:
#  -0.895194  -0.960283  -0.987055  -1.0       -0.968447  -1.0       -0.992712  -1.0       -1.0       -0.964975
#   0.438019   0.504086   0.623425   0.671433   0.475873   0.551948   0.563323   0.556987   0.577228   0.553053
#  -0.409617  -0.422937  -0.586731  -0.638229  -0.317215  -0.373129  -0.468687  -0.368036  -0.464358  -0.48796
#   0.313867   0.336202   0.509322   0.56959    0.233927   0.291145   0.386832   0.286327   0.384007   0.402437
#  -0.257968  -0.286627  -0.482991  -0.557164  -0.184029  -0.241863  -0.3434    -0.236681  -0.341143  -0.359296
```

The result is a matrix of size `n×m`, where `n` is the number of cepstral coefficients and `m` is the number of time series.

## Reference

```bibtex
@INPROCEEDINGS{989529,
  author={Kalpakis, K. and Gada, D. and Puttagunta, V.},
  booktitle={Proceedings 2001 IEEE International Conference on Data Mining}, 
  title={Distance measures for effective clustering of ARIMA time-series}, 
  year={2001},
  volume={},
  number={},
  pages={273-280},
  doi={10.1109/ICDM.2001.989529}}
```

