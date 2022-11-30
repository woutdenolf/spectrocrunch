# -*- coding: utf-8 -*-

from peewee import *
import urllib
import tempfile
import os

from contextlib import contextmanager
from sshtunnel import SSHTunnelForwarder
import traceback


db = MySQLDatabase(
    "cod", host="127.0.0.1", user="cod_reader", port=3308, connect_timeout=10000
)

# Get
# ssh wout@axil1.ua.ac.be -L 3307:www.crystallography.net:3306 -N &
# python -m pwiz cod -e mysql -u cod_reader -H 127.0.0.1 -p 3307
#
# mysql -ucod_reader -h 127.0.0.1 -P 3307
#   SELECT DATABASE();
#   USE cod;
#   SHOW TABLES;
#   DESCRIBE data;


class BaseModel(Model):
    class Meta:
        database = db


class Data(BaseModel):
    rfsqd = FloatField(db_column="RFsqd", null=True)
    ri = FloatField(db_column="RI", null=True)
    rall = FloatField(db_column="Rall", null=True)
    robs = FloatField(db_column="Robs", null=True)
    rref = FloatField(db_column="Rref", null=True)
    z = IntegerField(db_column="Z", index=True, null=True)
    zprime = FloatField(db_column="Zprime", index=True, null=True)
    a = FloatField(index=True, null=True)
    acce_code = CharField(index=True, null=True)
    alpha = FloatField(index=True, null=True)
    authors = TextField(null=True)
    b = FloatField(index=True, null=True)
    beta = FloatField(index=True, null=True)
    c = FloatField(index=True, null=True)
    calcformula = CharField(index=True, null=True)
    cellformula = CharField(null=True)
    cellpressure = FloatField(null=True)
    celltemp = FloatField(null=True)
    chemname = CharField(index=True, null=True)
    commonname = CharField(index=True, null=True)
    compoundsource = CharField(null=True)
    date = DateField(index=True, null=True)
    diffrpressure = FloatField(null=True)
    diffrtemp = FloatField(null=True)
    doi = CharField(index=True, null=True)
    duplicateof = IntegerField(null=True)
    file = PrimaryKeyField()
    firstpage = CharField(null=True)
    flags = CharField(null=True)
    formula = CharField(index=True, null=True)
    gamma = FloatField(index=True, null=True)
    gofall = FloatField(null=True)
    gofgt = FloatField(null=True)
    gofobs = FloatField(null=True)
    issue = CharField(null=True)
    journal = CharField(index=True, null=True)
    lastpage = CharField(null=True)
    method = CharField(index=True, null=True)
    mineral = CharField(index=True, null=True)
    nel = CharField(index=True, null=True)
    onhold = DateField(null=True)
    optimal = IntegerField(null=True)
    pressurehist = CharField(null=True)
    radsymbol = CharField(db_column="radSymbol", null=True)
    radtype = CharField(db_column="radType", null=True)
    radiation = CharField(null=True)
    sg = CharField(index=True, null=True)
    sghall = CharField(db_column="sgHall", index=True, null=True)
    siga = FloatField(null=True)
    sigalpha = FloatField(null=True)
    sigb = FloatField(null=True)
    sigbeta = FloatField(null=True)
    sigc = FloatField(null=True)
    sigcellpressure = FloatField(null=True)
    sigcelltemp = FloatField(null=True)
    sigdiffrpressure = FloatField(null=True)
    sigdiffrtemp = FloatField(null=True)
    siggamma = FloatField(null=True)
    sigvol = FloatField(null=True)
    status = CharField(null=True)
    svnrevision = IntegerField(index=True, null=True)
    text = TextField(index=True)
    thermalhist = CharField(null=True)
    time = TimeField(index=True, null=True)
    title = TextField(null=True)
    vol = FloatField(index=True, null=True)
    volume = IntegerField(null=True)
    wrall = FloatField(db_column="wRall", null=True)
    wrobs = FloatField(db_column="wRobs", null=True)
    wrref = FloatField(db_column="wRref", null=True)
    wavelength = FloatField(null=True)
    year = IntegerField(null=True)

    class Meta:
        db_table = "data"
        indexes = ((("mineral", "chemname", "commonname"), False),)

    def __str__(self):
        ret = "{} ({})\n".format(self.mineral, self.commonname)
        ret += "{} ({})\n".format(self.formula, self.chemname)
        ret += "{} ({} {} {} {} {} {})\n".format(
            self.sg, self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
        ret += "P = {} kPa, T = {} K\n".format(self.diffrpressure, self.diffrtemp)
        ret += "P = {} kPa, T = {} K\n".format(self.cellpressure, self.celltemp)
        ret += "{} ({})\n".format(self.authors, self.year)
        ret += "https://doi.org/{}\n".format(self.doi)
        return ret

    @staticmethod
    def sap(p):
        if p is None:
            return True
        a = 0.9  # atm
        b = 1.1
        a *= 101.325  # kPa
        b *= 101.325
        return p >= a and p <= b

    @staticmethod
    def sat(t):
        if t is None:
            return True
        a = 15  # celcius
        b = 30
        a += 273.15  # kelvin
        b += 273.15
        return t >= a and t <= b

    def satp(self):
        return (
            self.sap(self.diffrpressure)
            and self.sap(self.cellpressure)
            and self.sat(self.diffrtemp)
            and self.sat(self.celltemp)
        )

    @property
    def filename(self):
        return os.path.join("{}.cif".format(self.file))

    @property
    def path(self):
        return os.path.join(tempfile.gettempdir(), "spectrocrunch", "cif")

    @property
    def resourcename(self):
        return os.path.join(self.path, self.filename)

    @property
    def url(self):
        return "http://www.crystallography.net/cod/{}.cif".format(self.file)

    def download(self):
        filename = self.resourcename
        if not os.path.isfile(filename):
            path = self.path
            if not os.path.exists(path):
                os.makedirs(path)
            ciffile = urllib.URLopener()
            ciffile.retrieve(self.url, filename)

    @classmethod
    def namequery(cls, name):
        return (
            cls.select()
            .where(
                cls.mineral == name or cls.commonname == name or cls.chemname == name
            )
            .order_by(cls.year.desc())
        )


@contextmanager
def codtunnel():
    server = SSHTunnelForwarder(
        ssh_address_or_host=("axil1.ua.ac.be", 22),
        ssh_username="wout",
        ssh_pkey="/users/denolf/.ssh/id_rsa",
        remote_bind_address=("www.crystallography.net", 3306),
        local_bind_address=("127.0.0.1", 3308),
    )
    try:
        server.start()
        yield
    except:
        print(traceback.format_exc())

    server.stop()


if __name__ == "__main__":

    with codtunnel():

        query = Data.namequery("copper acetate")
        # for entry in query:
        #    print(entry)
        for entry in query:
            if entry.satp():
                print(entry)
                entry.download()
                break
