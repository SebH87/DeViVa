import pandas as pd

inf = pd.read_csv("datasub.txt", sep="\t",index_col=0, usecols=[0, 5, 8])
inf["ANN.GENE"] = inf["ANN.GENE"].fillna("?")
inf["ANN.AA"] = inf["ANN.AA"].fillna("?")
inf["AA_mut"] = inf["ANN.GENE"] + ":" + inf["ANN.AA"].apply(lambda x: str(x).strip("p."))
inf["label"] = inf.index + " (" + inf["AA_mut"] + ")"
inf["label"] = inf["label"].apply(lambda x: x.replace(" (?:?)", ""))
inf.to_csv("aa_label.txt", sep="\t")
