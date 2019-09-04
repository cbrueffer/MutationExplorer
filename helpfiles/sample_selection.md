# Sample Selection

Filter the overall sample list that survival analysis is conducted for. Filtering can be formed based on the treatment patients received and the status of several clinical biomarkers.  

### Treatment Group

Select a subgroup of patient samples treated with a specific combination three treatment regimens: endocrine therapy (Endo, e.g. tamoxifen), chemo therapy (Chemo), and anti HER2 therapy (HER2, e.g. trastuzumab). Available treatment combinations:

- Any
- Endocrine + any
- Chemo + any
- HER2 + any
- Endocrine only
- Chemo only
- HER2 only
- Endo + Chemo + any
- Endo + HER2 + any
- Cyto + HER2 + any
- Endo + Chemo only
- Endo + HER2 only
- Cyto + HER2 only
- Endo + Chemo + HER2
- No Systemic Treatment


### Biomarkers

| Biomarker            | Description                                                                                      |
|:-------------------- |:------------------------------------------------------------------------------------------------ |
| Histological Type    | Histological Type: Ductal, Lobular, or Other (mixed ductal/lobular, medullary, etc)              |
| ER                   | Estrogen receptor: Positive for &GreaterEqual;1%/10% IHC staining cells, otherwise Negative)     |
| PgR                  | Progesterone receptor: Positive for &GreaterEqual;1%/10% IHC staining cells, otherwise Negative) |
| HER2                 | Epidermal growth factor receptor 2 (HER2/ERBB2)                                                  |
| Ki67                 | Ki67 progression marker (staining cutoffs varying, depending on treating hospital)               |
| NHG                  | Nottingham histological grade                                                                    |
| PAM50                | Intrinsic breast cancer subtype based on the PAM50 gene list                                     |

The cutoffs ER and PgR can be either 1% (international standard) or 10% (Swedish standard).
