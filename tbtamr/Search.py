import pandas as pd
import pathlib,json
from tabulate import tabulate
from unidecode import unidecode

from CustomLog import logger

def check_file(pth) -> bool:

    if pathlib.Path(pth).exists():
        return True
    else:
        logger.critical(f"{pth} does not exist. Please try again.")
        raise SystemExit

def load_catalogue(pth) -> pd.DataFrame:

    if check_file(pth = pth):
        tab = pd.read_csv(pth, sep = None, engine = "python")
        return tab

def load_cfg(pth):

    if check_file(pth = pth):
        with open(f"{pth}", "r") as j:
            cfg = json.load(j)
            return cfg

def query_catalogue(tab, cfg, query):
    _results = []
    for q in query:
        r = tab[tab[cfg['variant_col']].str.contains(q)]
        # print(r)
        if not r.empty:
            _results.append(r)
    if _results != []:
        res = pd.concat(_results)
        vals = res.shape[0]
        cols = list(res.columns)
        _d = {key: [unidecode(key)] for key in cols}
        for row in res.iterrows():
            for col in cols:
                _d[col].append(row[1][col])
        res_list = [_d[col] for col in cols]
        mw = [30]
        for l in range(vals):
            mw.append(15)
        print(tabulate(res_list, tablefmt="simple_grid", maxcolwidths=mw))
    else:
        logger.warning(f"No results found for your query : {'|'.join(query)}")

def search(
            config,
            catalogue,
            query
            ):

    tab = load_catalogue(pth = catalogue)
    cfg = load_cfg(pth = config)
    query_catalogue(tab=tab, cfg=cfg, query=query)
