# p-B

These files are the R code and plot data files used for calculating species population turnover of biomass and abundance in a tropical forest plot.

**Paper**: Takkashi S. Kohyama, Matthew D. Potts, Tetsuo I. Kohyama, Kaoru Niiyama, Tze Leong Yao, Stuart J. Davies, & Douglas Sheil (2020) Trade-off between standing biomass and productivity in species-rich tropical forest: evidence, explanations and implications. *Journal of Ecology* in press.

## Contents

This dataset is a processed subset of the original dataset used in our analysis of biomass turnover across tree populations as demonstrated in the main paper. Readers interested in using the Pasoh 50-ha plot data for purposes other than reviewing our analysis are advised to contact the [Forest Research Institute Malaysia (FRIM)](https://www.frim.gov.my) and the [Center for Tropical Forest Science-Forest Global Earth Observatory (CTFS-Forest GEO)](https://forestgeo.si.edu), Smithsonian Tropical Research Institute.

* **turnover_pfr.r** - R code for calculating species-specific structural data and turnover rates.

* **data/pfr_observed.csv.gz** - Pasoh 50-ha plot data, Peninsular Malaysia, for ca. 1990 and ca. 2000 censuses: Cd for species code, (x, y) for coordinates (m) of stem location, dbh1 and dbh2 (cm) for stem diameter in 1990 and 2000, respectively, and t (year) for inter-census duration.

* **data/pfr_identity_free.csv.gz** - Pasoh 50-ha plot data (as pfr_observed.csv), where species identity and location are replaced between two stems with closest dbh1.

## License

The data files are licensed under the [Creative Commons Attribution 4.0 International License](LICENSE-CC-BY), and the R code is licensed under the [MIT License](LICENSE-MIT).
