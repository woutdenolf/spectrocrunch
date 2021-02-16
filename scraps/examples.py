# -*- coding: utf-8 -*-

if __name__ == "__main__":
    import glob
    from PyMca5.PyMcaIO import EdfFile

    files = glob.glob(
        "/data/id21/inhouse/15apr/Hiram/results/fXANES5/fXANES5_data/*.edf"
    )
    for f in files:
        h = EdfFile.EdfFile(f)
        data = h.GetData(0)
        print(f, data.shape)
