# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from ..utils.Enum import Enum
from ..utils import instance

import ast
import numpy as np
import pandas as pd
import scipy.stats
import xlsxwriter
import xlsxwriter.utility as xlsxutils
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__) 

class Writer(object):

    def __init__(self,filename=None):
        self.filename = filename
        self.open_file = None
        
    def __enter__(self):
        if self.filename:
            self.open_file = pd.ExcelWriter(self.filename,engine='xlsxwriter')
        else:
            self.open_file = None
        return self
        
    def __exit__(self,exc_type, exc_val, exc_tb):
        if self.open_file:
            self.open_file.save()
        self.open_file = None
        
    @property
    def handle(self):
        return self.open_file

    @property
    def workbook(self):
        return self.handle.book


class DataFrame(object):
    indextypes = Enum(['index','name','xlsindex','xlscell'])
    cellformats = Enum(['good','bad','select','scientific'])
    critereatypes = Enum(['between','notbetween','condition','colorscale'])
    priorities = Enum(['row','column'])
    
    def __init__(self,writer=None,sheet_name="Sheet1",priority="row",df=None,rowlevels=None,columnlevels=None,**kwargs):
        self.writer = writer
        self.sheet_name = sheet_name
        self.priority = self.priorities(priority)
        
        if rowlevels:
            kwargs["index"] = pd.MultiIndex.from_tuples([], names=rowlevels)
        if columnlevels:
            kwargs["index"] = pd.MultiIndex.from_tuples([], names=columnlevels)
        
        if df is None:
            self.df = pd.DataFrame(**kwargs)
        else:
            self.df = df
            
        self.criteria_column = OrderedDict()
        self.criteria_row = OrderedDict()
        self.formats_column = OrderedDict()
        self.formats_row = OrderedDict()
        self.formulae_column = OrderedDict()
        self.formulae_row = OrderedDict()
        self.formulae_cell = OrderedDict()
        self.headerfmt = {"bold":True,"bg_color":'#ffffff',"font_color":"#000000"}
    
    @classmethod
    def fromexcel(cls,filename,sheet_name=None):
        data = pd.read_excel(filename,sheet_name=sheet_name)
        if sheet_name is None:
            return [cls(df=df,sheet_name=sheet_name) for sheet_name,df in data.items()]
        else:
            return [cls(df=data,sheet_name=sheet_name)]
    
    def addvalue(self,row,column,data):
        self._remove_formulae(row=row,column=column)
        self._addvalue(row,column,data)
        self._reapply_formulae(row=row,column=column)
 
    def addrow(self,row,data):
        self._remove_formulae(row=row)
        self._addrow(row,data)
        self._reapply_formulae(row=row)
         
    def addcolumn(self,column,data):
        self._remove_formulae(column=column)
        self._addcolumn(column,data)
        self._reapply_formulae(column=column)

    def concat(self,df):
        self.df = pd.concat([self.df,df],ignore_index=False,sort=False)

    def _addvalue(self,row,column,data):
        self.df.at[row,column] = data
        
    def _addrow(self,row,data):
        if isinstance(data,dict):
            for column,value in data.items():
                self.df.at[row,column] = value
        else:
            self.df.at[row,:] = data

    def _addcolumn(self,column,data):
        if isinstance(data,dict):
            for rowindew,value in data.items():
                self.df.at[rowindew,column] = value
        else:
            self.df.at[:,column] = data

    def _reapply_formulae(self,row=None,column=None):
        self.apply_formulae()
        # Remark: differentiated reapply requires testing recursive formulae which is time consuming
                    
    def _remove_formulae(self,row=None,column=None):
        if row is not None and column is not None:
            pass
        elif row is not None:
            self.formulae_cell = {(r,c):v for (r,c),v in self.formulae_cell.items() if r!=row}
        elif column is not None:
            self.formulae_cell = {(r,c):v for (r,c),v in self.formulae_cell.items() if c!=column}

    def addcolumn_formula(self,column,formula,columns):
        if not instance.isarray(columns):
            columns = [columns]  
        if column in columns:
            rhside = formula.format(*columns)
            raise RuntimeError("Self referencing formula: {} = {}".format(column,rhside))
        self.formulae_column[column] = formula,columns
        self._apply_column_formula(column,formula,columns)
        self._reapply_formulae(column=column)
        
    def addrow_formula(self,row,formula,rows):
        if not instance.isarray(rows):
            rows = [rows]
        if row in rows:
            rhside = formula.format(*rows)
            raise RuntimeError("Self referencing formula: {} = {}".format(row,rhside))
        self.formulae_row[row] = formula,rows
        self._apply_row_formula(row,formula,rows)
        self._reapply_formulae(row=row)
        
    def addcell_formula(self,row,column,formula,rows,columns):
        if not instance.isarray(rows):
            rows = [rows]
        if not instance.isarray(columns):
            columns = [columns]  
        if row in rows or column in columns:
            rhside = formula.format(*[rc for rc in zip(rows,columns)])
            raise RuntimeError("Self referencing formula: {} = {}".format(row,column,rhside))
        self.formulae_cell[(row,column)] = formula,rows,columns
        self._apply_cell_formula(row,column,formula,rows,columns)
        self._reapply_formulae(row=row,column=column)
        
    def _apply_cell_formula(self,row,column,formula,rows,columns):
        data = self._formula_eval(formula,rows,columns,column=False)
        self._addvalue(row,column,data)

    def _apply_column_formula(self,column,formula,columns,indices=None):
        if indices is None:
            data = self._formula_eval_range(formula,columns,column=True)
            self._addcolumn(column,data)
        else:
            for row in instance.arrayit(indices):
                data = self._formula_eval(formula,columns,row,column=True)
                self._addvalue(row,column,data)
        
    def _apply_row_formula(self,row,formula,indices,columns=None):
        if columns is None:
            data = self._formula_eval_range(formula,indices,column=False)
            self._addrow(row,data)
        else:
            for column in instance.arrayit(columns):
                data = self._formula_eval(formula,indices,column,column=False)
                self._addvalue(row,column,data)
            
    def _apply_column_formulae(self,indices=None):
        for column,args in self.formulae_column.items():
            self._apply_column_formula(column,*args,indices=indices)
    
    def _apply_row_formulae(self,columns=None):
        for row,args in self.formulae_row.items():
            self._apply_row_formula(row,*args,columns=columns)

    def _apply_cell_formulae(self,row=None,column=None):
        if row is not None and column is not None:
            fapply = lambda r,c: r==row and c==column
        elif row is not None:
            fapply = lambda r,c: r==row
        elif column is not None:
            fapply = lambda r,c: c==column
        else:
            fapply = lambda r,c: True
            
        for (row,column),args in self.formulae_cell.items():
            if fapply(row,column):
                self._apply_cell_formula(row,column,*args)

    def apply_formulae(self):
        if self.priority==self.priorities.row:
            self._apply_column_formulae()
            self._apply_row_formulae()
        else:
            self._apply_row_formulae()
            self._apply_column_formulae()
        self._apply_cell_formulae()
        
    @classmethod
    def _formula_argfunc(cls,x,y,swap=False):
        if swap:
            x,y = y,x
        if instance.isstring(x):
            x = "\"{}\"".format(x)
        if instance.isstring(y):
            y = "\"{}\"".format(y)
        return "self.df.at[{},{}]".format(x,y) # row,column
        
    @classmethod
    def _formula_argfunc_range(cls,x,column=True):
        if column:
            add = ""
        else:
            add = ".loc"
        if instance.isstring(x):
            return "self.df{}[\"{}\"]".format(add,x)
        else:
            return "self.df{}[{}]".format(add,x)
    
    def _eval_selfctx(self,expr):
        expr = ast.parse(expr, mode='eval')
        code = compile(expr, filename="<ast>", mode="eval")
        return eval(code, globals(), locals())
        #return eval(expr, globals(), locals())
        
    def _formula_eval_range(self,expr,args,column=True):
        # column refers to args
        expr = expr.format(*[self._formula_argfunc_range(arg,column=column) for arg in args])
        return self._eval_selfctx(expr)
        
    def _formula_eval(self,expr,args,others,column=True):
        # column refers to args
        if not instance.isarray(args):
            args = [args]
        if not instance.isarray(others):
            others = [others]*len(args)
        expr = expr.format(*[self._formula_argfunc(arg,other,swap=column) for arg,other in zip(args,others)])
        return self._eval_selfctx(expr)

    def addrowformat(self,index,fmt=None):
        self._addformat(index,fmt=fmt,column=False)
    
    def addcolumnformat(self,index,fmt=None):
        self._addformat(index,fmt=fmt,column=True)
        
    def _addformat(self,index,fmt=None,column=True):
        if fmt:
            if column:
                if index in self.formats_column:
                    fmt = self._append_format(self.formats_column[index]["fmt"],fmt)
                self.formats_column[index] = {"fmt":fmt}
            else:
                if index in self.formats_row:
                    fmt = self._append_format(self.formats_row[index]["fmt"],fmt)
                self.formats_row[index] = {"fmt":fmt}
                
    def _append_format(self,fmt,add):
        if not instance.isarray(fmt):
            fmt = [fmt]
        if instance.isarray(add):
            fmt.extend(add)
        else:
            fmt.append(add)
        return fmt
             
    def addcriterium(self,index,crit=None,column=True):
        if column:
            self.criteria_column[index] = crit
        else:
            self.criteria_row[index] = crit
            
    def addcriterium_outliers(self,index,alpha=0.9,column=True,out=True,sided="double"):
        crit = {}
        if sided=="single":
            m = scipy.stats.norm.ppf(alpha)
            if out:
                crit["mode"] = self.critereatypes.condition
                crit["operator"] = ">"
            else:
                crit["mode"] = self.critereatypes.condition
                crit["operator"] = "<="
        else:
            m = np.array(scipy.stats.norm.interval(alpha))
            if out:
                crit["mode"] = self.critereatypes.notbetween
            else:
                crit["mode"] = self.critereatypes.between
                
        crit["func"] = lambda x: np.nanmean(x)+m*np.nanstd(x)
        
        if out:
            crit["fmt"] = self.cellformats.bad
        else:
            crit["fmt"] = self.cellformats.good

        self.addcriterium(index,crit,column=column)

    def addcriterium_colorscale(self,index,fmt=None,column=True,colors=None):
        if colors is None:
            #green,yellow,red
            colors = ['#63be7b','#ffeb84','#f8696b']
            
        if len(colors)==3:
            cond = {'type': '3_color_scale',
                    'min_color': colors[0],
                    'mid_color': colors[1],
                    'max_color': colors[2]}
        else:
            cond = {'type': '2_color_scale',
                    'min_color': colors[0],
                    'max_color': colors[-1]}
  
        crit = {"mode":self.critereatypes.colorscale,"cond":cond}
        if fmt:
            crit["fmt"] = fmt
        self.addcriterium(index,crit,column=column)
    
    @property
    def handle(self):
        if self.writer:
            return self.writer.handle
        else:
            return None
        
    @property
    def workbook(self):
        return self.writer.workbook

    @property
    def worksheet(self):
        return self.handle.sheets[self.sheet_name]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.save()
    
    def save(self):
        writer = self.handle
        if writer and not self.df.empty:
            self.df.to_excel(writer,sheet_name=self.sheet_name)
            self._xls_save_formats()
            self._xls_save_colwidths()

    def sortrows(self,level=0):
        self.df.sort_index(axis=0,level=level,sort_remaining=True,inplace=True)

    def sortcolumns(self,level=0):
        self.df.sort_index(axis=1,level=level,sort_remaining=True,inplace=True)
     
    def sort(self,rowlevel=0,columnlevel=0):
        self.sortrows(level=rowlevel)
        self.sortcolumns(level=columnlevel)
    
    def _xls_save_colwidths(self):
        worksheet = self.worksheet
        
        r0,c0 = self._xls_rowoff,self._xls_coloff
        
        # Width of index columns:
        width = 15
        
        for c in range(c0):
            if c0==1:
                width = max(width,max(len(rowindex) for rowindex in self.df.index))
            else:
                width = max(width,max(len(rowindex[c]) for rowindex in self.df.index))
            hname = self.df.index.names[c]
            if hname:
                width = max(width,len(hname))
                        
            worksheet.set_column(c,c,width)

        # Width of data columns:
        for i,colindex in enumerate(self.df.columns,c0):
            if r0==1:
                width = max(len(colindex),15)
            else:
                width = max(max(len(c) for c in colindex),15)
            worksheet.set_column(i,i,width)
        
        # Freeze index and column names:
        worksheet.freeze_panes(r0,c0)

    def _xls_save_formats(self):
        for index,crit in self.criteria_column.items():
            self._xls_apply_criterium(index,crit,column=True)
        for index,crit in self.criteria_row.items():
            self._xls_apply_criterium(index,crit,column=False)
        for index,fmt in self.formats_column.items():
            self._xls_apply_format(index,fmt,column=True)
        for index,fmt in self.formats_row.items():
            self._xls_apply_format(index,fmt,column=False)
        self._xls_save_headerformats()
        
    def _xls_save_headerformats(self):
        fmt = self._xls_add_format(self.headerfmt)
        for c in range(self._xls_coloff):
            self.worksheet.set_column(c,c,None,fmt)
        for r in range(self._xls_rowoff):
            self.worksheet.set_row(r,None,fmt)
        
    def _xls_add_format(self,fmts):
        fmtdict = self._xls_parse_format(fmts)
        if fmtdict:
            return self.workbook.add_format(fmtdict)
        else:
            return None
            
    def _xls_parse_format(self,fmts):
        fmtdict = {}
        if not instance.isarray(fmts):
            fmts = [fmts]
        for fmt in fmts:
            if fmt==self.cellformats.good:
                fmtdict["bg_color"] = '#C6EFCE'
                fmtdict["font_color"] = '#006100'
            elif fmt==self.cellformats.bad:
                fmtdict["bg_color"] = '#FFC7CE'
                fmtdict["font_color"] = '#9C0006'
            elif fmt==self.cellformats.select:
                fmtdict["bg_color"] = '#FFEB9C'
                fmtdict["font_color"] = '#9C6500'
            elif fmt==self.cellformats.scientific:
                fmtdict["num_format"] = "0.00E+00"
            else:
                if isinstance(fmt,dict):
                    fmtdict.update(fmt)
                elif fmt:
                    fmtdict["bg_color"] = '#FFFFFF'
                    fmtdict["font_color"] = '#000000'
        return fmtdict

    def _xls_apply_criterium(self,index,crit,column=True):
        fromtype = crit.get("fromtype",self.indextypes.name)
        mode = crit.get("mode",self.critereatypes.between)
        fmt = crit.get("fmt",self.cellformats.select)

        # Select data and excel range
        if column:
            head,ran = self._xls_column(index,fromtype)
            name = self._convert_colindex(index,fromtype,self.indextypes.name)
            values = self.df[name].values
        else:
            head,ran = self._xls_row(index,fromtype)
            name = self._convert_rowindex(index,fromtype,self.indextypes.name)
            values = self.df.loc[name].values
            
        # Add condition
        if mode==self.critereatypes.between or mode==self.critereatypes.notbetween:
            if mode==self.critereatypes.between:
                omi,oma,criteria = "<",">","between"
            else:
                omi,oma,criteria = ">","<","not between"
                
            func = crit.get("func",lambda x:(np.nanmin(x),np.nanmax(x)))
            mi,ma = func(values)
            if mi == -np.inf:
                cond = {"type":"cell","criteria":omi,"value":ma}
            elif ma == np.inf:
                cond = {"type":"cell","criteria":oma,"value":mi}
            else:
                cond = {"type":"cell","criteria":criteria,"minimum":mi,"maximum":ma}
        elif mode==self.critereatypes.condition:
            func = crit.get("func",lambda x:np.nanmean(x))
            op = crit.get("operator","==")
            cond = {"type":"cell","criteria":op,"value":func(values)}
        elif mode==self.critereatypes.colorscale:
            cond = crit.get("cond",{'type': '3_color_scale'})
        else:
            cond = {"type":"cell","criteria":">","value":-1e9}
        
        cond["format"] = self._xls_add_format(fmt)
        self.worksheet.conditional_format(ran, cond)

    def _xls_apply_format(self,index,params,column=True):
        fromtype = params.get("fromtype",self.indextypes.name)
        fmt = params.get("fmt",None)

        fmt = self._xls_add_format(fmt)
        if column:
            coli = self._convert_colindex(index,fromtype,self.indextypes.xlsindex)
            self.worksheet.set_column(coli,coli,None,fmt)
        else:
            rowi = self._convert_rowindex(index,fromtype,self.indextypes.xlsindex)
            self.worksheet.set_row(rowi,None,fmt)

    @property
    def _xls_coloff(self):
        return len(self.df.index.names)

    @property
    def _xls_rowoff(self):
        return len(self.df.columns.names)

    @property
    def _xls_ncol(self):
        return self.df.shape[1]

    @property
    def _xls_nrow(self):
        return self.df.shape[0]

    def _convert_singleindex(self,x,fromtype,totype,row=True):
        if fromtype==totype:
            return x
        
        if fromtype==self.indextypes.index:
            if totype==self.indextypes.name:
                if row:
                    y = self.df.index[x]
                else:
                    y = self.df.columns[x]
            else: # xlsindex, xlscell
                if row:
                    y = x+self._xls_rowoff
                else:
                    y = x+self._xls_coloff
        elif fromtype==self.indextypes.name:
            if row:
                y = self.df.index.get_loc(x)
            else:
                y = self.df.columns.get_loc(x)
            y = self._convert_singleindex(y,self.indextypes.index,totype,row=row)
        else: # xlsindex, xlscell
            if row:
                y = x-self._xls_rowoff
            else:
                y = x-self._xls_coloff
            y = self._convert_singleindex(y,indextypes.index,totype,row=row)

        return y

    def _convert_rowindex(self,*args):
        return self._convert_singleindex(*args,row=True)

    def _convert_colindex(self,*args):
        return self._convert_singleindex(*args,row=False)

    def _convert_cell(self,row,col,fromtype,totype,rowabs=False,colabs=False):
        if fromtype==totype:
            return row,col
            
        if fromtype==self.indextypes.xlscell:
            row,col = xlsxutils.xl_cell_to_rowcol(row,col,rowabs,colabs)
            fromtype = self.indextypes.xlsindex
        
        row = self._convert_rowindex(row,fromtype,totype)
        col = self._convert_colindex(col,fromtype,totype)
        
        if totype==self.indextypes.xlscell:
            row,col = xlsxutils.xl_rowcol_to_cell(row,col,rowabs,colabs)

        return row,col

    def _xls_column(self,col,fromtype,absolute=True):
        coli = self._convert_colindex(col,fromtype,self.indextypes.xlscell)

        head = xlsxutils.xl_rowcol_to_cell(0,coli,True,True)
        if absolute:
            ran = xlsxutils.xl_range_abs(self._xls_rowoff,coli,self._xls_nrow+self._xls_rowoff-1,coli)
        else:
            ran = xlsxutils.xl_range(self._xls_rowoff,coli,self._xls_nrow+self._xls_rowoff-1,coli)
        return head,ran

    def _xls_row(self,row,fromtype,absolute=True):
        rowi = self._convert_rowindex(row,fromtype,self.indextypes.xlscell)

        head = xlsxutils.xl_rowcol_to_cell(rowi,0,True,True)
        if absolute:
            ran = xlsxutils.xl_range_abs(rowi,self._xls_coloff,rowi,self._xls_ncol+self._xls_coloff-1)
        else:
            ran = xlsxutils.xl_range(rowi,self._xls_coloff,rowi,self._xls_ncol+self._xls_coloff-1)
        return head,ran


