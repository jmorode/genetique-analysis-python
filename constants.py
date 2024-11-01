import numpy as np

# ML Relate
reference_array = np.array(
    [
        ["-", np.nan, np.nan, np.nan, np.nan, np.nan],
        ["U", "-", np.nan, np.nan, np.nan, np.nan],
        ["U", "U", "-", np.nan, np.nan, np.nan],
        ["PO", "PO", "U", "-", np.nan, np.nan],
        ["PO", "PO", "U", "FS", "-", np.nan],
        ["U", "PO", "PO", "HS", "HS", "-"],
    ],
    dtype=object,
)

order_inds_family = {"P1": 1, "M": 2, "P2": 3, "F1": 4, "F2": 5, "DF": 6}

nb_ind_relation = {"PO": 6, "U": 6, "FS": 1, "HS": 2}


# Generation
rename_parents = dict({"Ind1": "P1", "Ind2": "M", "Ind3": "P2"})
