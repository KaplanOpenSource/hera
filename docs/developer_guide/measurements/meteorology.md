# Meteorology toolkits

## lowFreqToolKit

The toolkit follows the standard analysis + presentation pattern:

| Layer | Class | Key methods |
|-------|-------|------------|
| Toolkit | `lowFreqToolKit` | `getDataSourceData`, standard datasource API |
| Analysis | `lowFreqAnalysis` | `addDatesColumns`, `calcHourlyDist`, `resampleSecondMoments`, `matchDataWithOther` |
| Presentation | `lowFreqPresentation` | `dailyPlots.plotScatter`, `dailyPlots.plotDaily`, `seasonalPlots.plotProbContourf_bySeason`, `seasonalPlots.plotSeasonalHourly` |

## HighFreqToolKit

High-frequency data uses a calculator-based analysis architecture:

| Layer | Class | Purpose |
|-------|-------|---------|
| Toolkit | `HighFreqToolKit` | Data source management, station metadata |
| Analysis | `RawdataAnalysis` | Dispatches to specialized calculators |
| Calculator | `AbstractCalculator` | Base class for analysis calculators |
| Calculator | `MeanDataCalculator` | Mean statistics (wind speed, direction, temperature) |
| Calculator | `TurbulenceStatistics` | Second-order moments, eddy covariance |
| Parsers | `CampbellBinary`, `TOA5` | Raw data format parsers |
