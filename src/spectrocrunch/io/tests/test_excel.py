import unittest
import numpy as np
import pandas as pd
import random
from testfixtures import TempDirectory
import os

from .. import excel


class test_excel(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()
        self.nrow, self.ncol = 8, 6
        self._rows = [chr(c + ord("a")) for c in range(0, self.nrow)]
        self._columns = [chr(c + ord("a")) for c in range(13, 13 + self.ncol)]

    def rows(self, n=None):
        rows = [chr(c + ord("a")) for c in range(0, self.nrow)]
        if n:
            np.random.shuffle(rows)
            rows = rows[:n]
        return rows

    def columns(self, n=None):
        columns = [chr(c + ord("a")) for c in range(13, 13 + self.ncol)]
        if n:
            np.random.shuffle(columns)
            columns = columns[:n]
        return columns

    def tearDown(self):
        self.dir.cleanup()

    def _generate_df(self, writer=None):
        data = np.random.random(self.nrow * self.ncol).reshape((self.nrow, self.ncol))
        priority = random.choice(list(excel.DataFrame.priorities))
        df = excel.DataFrame(
            writer=writer,
            data=data,
            columns=self.columns(),
            index=self.rows(),
            priority=priority,
        )
        for i in range(5):
            m = random.choice([0, 1, 2])
            rows = self.rows(4)
            columns = self.columns(4)
            if m == 2:
                df.addrow_formula(rows[0], "({}+{}+{})/3.", rows[1:])
            elif m == 1:
                df.addcolumn_formula(columns[0], "({}+{}+{})/3.", columns[1:])
            else:
                df.addcell_formula(
                    rows[0], columns[0], "({}+{}+{})/3.", rows[1:], columns[1:]
                )

        return df

    def test_save(self):
        filename = os.path.join(self.dir.path, "test.xlsx")
        with excel.Writer(filename) as writer:
            df1 = self._generate_df(writer=writer)
            df1.save()
        df2 = next(iter(excel.DataFrame.fromexcel(filename, index_col=0).values()))
        np.testing.assert_array_almost_equal(df1.df.values, df2.df.values)

    def test_dataframe(self):
        data = np.arange(self.nrow * self.ncol).reshape((self.nrow, self.ncol))

        rows = self.rows()
        columns = self.columns()

        for equal in [True, False]:
            df2 = pd.DataFrame(data=data, columns=columns, index=rows)

            if equal:
                priority = excel.DataFrame.priorities.column
            else:
                priority = excel.DataFrame.priorities.row

            df1 = excel.DataFrame(columns=columns, index=rows, priority=priority)
            for index, row in zip(rows, data):
                df1.addrow(index, row)

            for i in [2, 3]:
                j, k = i - 2, i - 1
                df1.addrow_formula(rows[i], "({}+{})/2.", [rows[j], rows[k]])
                df2.loc[rows[i]] = (df2.loc[rows[j]] + df2.loc[rows[k]]) / 2

            for i in [2, 3]:
                j, k = i - 2, i - 1
                df1.addcolumn_formula(
                    columns[i], "({}+{})/3.", [columns[j], columns[k]]
                )
                df2[columns[i]] = (df2[columns[j]] + df2[columns[k]]) / 3

            i = 2
            df1.addcell_formula(
                rows[i], columns[i], "{}", rows[self.nrow - 1], columns[self.ncol - 1]
            )
            df2.loc[rows[i], columns[i]] = df2.loc[
                rows[self.nrow - 1], columns[self.ncol - 1]
            ]

            # dfequal = np.array_equal(df1.df.values,df2.values)
            # dfequal = df1.df.equals(df2)
            if equal:
                np.testing.assert_array_equal(df1.df.values, df2.values)
            else:
                np.testing.assert_raises(
                    AssertionError,
                    np.testing.assert_array_equal,
                    df1.df.values,
                    df2.values,
                )
