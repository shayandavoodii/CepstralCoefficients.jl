# Time Series Cepstral Clustering

This repo contains the Julia implementation of a study entitled "Distance Measures for Effective Clustering of ARIMA Time-Series"[[1](https://doi.org/10.1109/ICDM.2001.989529)]. The research proposes a distinguished approach to classify time series regarding their inner patterns. The authors used the Cepstral Analysis concept to find the patterns underneath the time series.

## How to use?

To use the implementation, one should pursue the following steps:

1. Install the package by running the following command in the Julia REPL:  

    ```julia
    using Pkg
    Pkg.add(url="https://github.com/shayandavoodii/TimeSeries-Cepstral-Clustering.git")
    ```

2. Import the package using the following command in the Julia REPL:

    ```julia
    using CepstralCoefficients
    ```

Afterward, one should use the [`cc`](https://github.com/shayandavoodii/TimeSeries-Cepstral-Clustering/blob/main/src/CepstralClustering.jl#L18-L44) function as the primary function of the implementation. The function above calculates the first `n` number of cepstral coefficients according to the given time series. The following methods for calculating the coefficients are available so far:

- Auto-Regressive Cepstral Coefficients
- Real Cepstral Coefficients

The former method can be employed by passing the `ARCepstral` type as the first argument to the `cc` function, and the latter can be utilized by passing the `RealCepstral` type instead.

## Example

A sequence of time series is required to test the implementation. In this case, the following assets are used: `["MSFT", "AAPL", "GOOG", "AMZN", "FB", "TSLA", "NVDA", "INTC", "CSCO", "ADBE"]`. In this example, the first 5 cepstral coefficients are calculated. The lag of $AR$ process is set to be `3`.

### Fetch data

```julia
using YFinance
tickers = ["MSFT", "AAPL", "GOOG", "AMZN", "META", "TSLA", "NVDA", "INTC", "CSCO", "ADBE"]
querry = [get_prices(ticker, startdt="2019-01-01", enddt="2020-01-01")["adjclose"] for ticker in tickers]
prices = stack(querry, dims=1);
```

Afterward, the [`cc`](https://github.com/shayandavoodii/TimeSeries-Cepstral-Clustering/blob/b586666d6764ac4e742cc07549c0247be30baa1b/src/CepstralClustering.jl#L12-L46) function is employed to calculate the cepstral coefficients.

### Calculate cepstral coefficients

```julia
n=5
p=3
model = ARCepstral(p)
cepscoefs = cc(model, prices, n)
# 5×10 Matrix{Float64}:
#  -0.895194  -0.960283  -0.987055  -1.0       -0.968447  -1.0       -0.992712  -1.0       -1.0       -0.964975
#   0.438019   0.504086   0.623425   0.671433   0.475873   0.551948   0.563323   0.556987   0.577228   0.553053
#  -0.409617  -0.422937  -0.586731  -0.638229  -0.317215  -0.373129  -0.468687  -0.368036  -0.464358  -0.48796
#   0.313867   0.336202   0.509322   0.56959    0.233927   0.291145   0.386832   0.286327   0.384007   0.402437
#  -0.257968  -0.286627  -0.482991  -0.557164  -0.184029  -0.241863  -0.3434    -0.236681  -0.341143  -0.359296
```

The result is a `n×m` matrix, where `n` is the number of cepstral coefficients and `m` is the number of time series. The next step is to perform clustering on the calculated cepstral coefficients.

### Clustering

In this regard, one can use the `clustering` function to perform Partition Around Medoids (PAM) clustering on the calculated cepstral coefficients. The function takes the cepstral coefficients as the first argument and the maximum number of clusters to be examined (in order to find the optimal number of clusters) as the second argument.

```julia
using Clustering
k = 4
clusters = cepsclustering(cepscoefs, k)
# 2-element Vector{Vector{Int64}}:
#  [1, 3, 5, 7, 9]
#  [2, 4, 6, 8, 10]
```

The result indicates that the 1st, 3rd, 5th, 7th, and 9th time series are in the first cluster and the rest are in the second cluster. In other words, the first cluster contains "MSFT", "GOOG", "META", "NVDA", and "CSCO" and the second cluster contains "AAPL", "AMZN", "FB", "TSLA", and "INTC". Note that, the result may vary each time due to the random nature of the PAM algorithm.

### Plotting

In order to probe the results, it is better to visualize it. In this subsection, the time series are plotted in two different color tones each of which represents a cluster.

```julia
colortones = ["#E46262" "#2C8EF6" "#CD6969" "#408BDA" "#B57070" "#5387BF" "#9E7777" "#6784A3" "#877E7E" "#7A8188"]
linestyle = [:solid :dash :solid :dash :solid :dash :solid :dash :solid :dash]
plot(
  prices',
  label=permutedims(tickers),
  title="Stock prices\nAnalysed data range: $(startdt) to $(enddt)",
  xlabel="Date",
  ylabel="Price",
  legend=:outerright,
  dpi=300,
  color=colortones,
  size=(1000, 400),
  linestyle=linestyle,
  bottom_margin=6mm,
  left_margin=6mm,
)
```
![img](https://github.com/shayandavoodii/TimeSeries-Cepstral-Clustering/blob/main/assets/StockPrices.png)

The results are not satisfactory, which is expected since the PAM clustering is inaccurate due to its random initialization. The random initialization may result in a nonoptimal solution. As seen in the figure above, the 'ABDE' and 'NVDA' series follow similar patterns but are in different clusters; this is surprising because the opposite was expected.

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
  doi={10.1109/ICDM.2001.989529}
}
```

