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

from ..common.Enum import Enum

import numpy as np
import pandas as pd
import scipy.stats
import xlsxwriter
import xlsxwriter.utility as xlsxutils

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
    cellformats = Enum(['good','bad','select'])
    critereatypes = Enum(['between','notbetween','condition','colorscale'])
    
    def __init__(self,writer=None,sheet_name="diagnostics"):
        self.writer = writer
        self.sheet_name = sheet_name
        self.df = pd.DataFrame()
        self.criteria = {}
        self.formats = {}
        
    def adddata(self,index,dic):
        for k,v in dic.items():
            self.df.at[index,k] = v

    def addcriterium_outliers(self,index,alpha=0.9,col=True,out=True,sided="double"):
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

        self.addcriterium(index,crit,col=col)

    def addcriterium_colorscale(self,index,col=True,colors=None):
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
        self.addcriterium(index,crit,col=col)
        
    def addcriterium(self,index,crit=None,col=True):
        crit["col"] = col
        self.criteria[index] = crit
    
    def addformat(self,index,fmt=None,col=True):
        self.formats[index] = {"fmt":fmt,"col":col}
    
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

    def _xls_save_colwidths(self):
        worksheet = self.worksheet
        width = max(len(row) for row in self.df.index)
        width = max(width,15)
        worksheet.set_column(0,0, width)
        for i,col in enumerate(self.df.columns):
            width = max(len(col),15)
            worksheet.set_column(i+1,i+1, width)
        worksheet.freeze_panes(1,1)

    def _xls_save_formats(self):
        for index,crit in self.criteria.items():
            self._xls_apply_criterium(index,crit)
        for index,fmt in self.formats.items():
            self._xls_apply_format(index,fmt)

    def _xls_add_format(self,fmt):
        return self.workbook.add_format(self._xls_parse_format(fmt))

    def _xls_parse_format(self,fmt):
        if isinstance(fmt,dict):
            return fmt
        else:
            name = fmt
            fmt = {}
            if name==self.cellformats.good:
                fmt["bg_color"] = '#C6EFCE'
                fmt["font_color"] = '#006100'
            elif name==self.cellformats.bad:
                fmt["bg_color"] = '#FFC7CE'
                fmt["font_color"] = '#9C0006'
            elif name==self.cellformats.select:
                fmt["bg_color"] = '#FFEB9C'
                fmt["font_color"] = '#9C6500'
            else:
                if not isinstance(name,dict):
                    fmt["bg_color"] = '#FFFFFF'
                    fmt["font_color"] = '#000000'
            return fmt

    def _xls_apply_criterium(self,index,crit):
        fromtype = crit.get("fromtype",self.indextypes.name)
        col = crit.get("col",True)
        mode = crit.get("mode",self.critereatypes.between)
        fmt = crit.get("fmt",self.cellformats.select)
        
        # Select data and excel range
        if col:
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

    def _xls_apply_format(self,index,params):
        fromtype = params.get("fromtype",self.indextypes.name)
        col = params.get("col",True)
        fmt = params.get("fmt",self.cellformats.none)

        fmt = self._xls_add_format(fmt)
        if col:
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
        return 1

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
                    
            y = self._convert_singleindex(y,self.indextypes.index,totype)
        else: # xlsindex, xlscell
            if row:
                y = x-self._xls_rowoff
            else:
                y = x-self._xls_coloff
            
            y = self._convert_singleindex(y,indextypes.index,totype)

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


