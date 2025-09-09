## Test environments
- Local: macOS Sequoia, R 4.5.1
- win-builder: devel / release / oldrelease (all OK)
- R-hub: windows, macOS-arm64, linux (all OK)

## R CMD check results
0 errors | 0 warnings | 0 notes
*(On local, `--as-cran` sometimes shows a NOTE: “CRAN incoming feasibility: New submission”. 
This is not a new submission but an update of the existing package.)*

## Notes for CRAN
- This is an **update release** of the existing CRAN package **darthtools**, not a new submission.  
- License remains MIT + file LICENSE.
- Updated code requires R (>= 4.1.0) due to use of base pipe (`|>`).  
