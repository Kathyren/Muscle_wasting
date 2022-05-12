import os

from eutils import Client
from eutils._internal.xmlfacades.pubmedarticleset import PubmedArticleSet


def _override_esearch(self, db, term, retmax):
    """query the esearch endpoint
    """
    esr = ESearchResult(self._qs.esearch({"db": db, "term": term, "retmax": retmax}))
    if esr.count > esr.retmax:
        print("NCBI found {esr.count} results, but we truncated the reply at {esr.retmax}"
              " results; see https://github.com/biocommons/eutils/issues/124/".format(esr=esr))
    return esr


import lxml.etree as le


def _elink(self, db, id):
    """query the elink endpoint
        """
    db = db.lower()
    xml = self._qs.elink({"db": db, "dbfrom": "pubmed", "linkname": "pubmed_pmc_refs", "id": str(id)})
    import xml.etree.ElementTree as ET
    xml = xml.decode('latin')
    root = ET.fromstring(xml)
    return root


Client.esearch = _override_esearch
Client.elink = _elink
import enum

from eutils._internal.xmlfacades.esearchresult import ESearchResult


class NCBIDatabases(enum.Enum):
    """
    This will store the name of the databases and
    their respective codes and stuff
    """

    def __new__(cls, *args, **kwds):
        value = len(cls.__members__) + 1
        obj = object.__new__(cls)
        obj._value_ = value
        return obj

    def __init__(self, database, main_attribute):
        self.database = database
        self.main_attribute = main_attribute

    Nucleotides = 'nuccore', 'gbseqs'
    Genes = 'gene', ''
    Pubmed = 'pubmed', '_xml_root'


class EutilsConnection():
    ec = None
    database = None

    def __init__(self, database: NCBIDatabases):
        self.database = database
        self.db = database.database
        self.table = database.main_attribute
        self.ec = Client(api_key=os.environ.get("NCBI_API_KEY", None))

    def fetch_queries_ids(self, term: str, db=None, retmax=100000) -> list:
        """
        This function will call the NCBI database "db" and look for the query specified on
        'term'. It could be as simple as "XM_537211". It will retrieve the query information
        with all the ids that match the search. If the get_first is True, then it will only get
        the first result.
        :param db:
        :param term:
        :param get_first:
        :return: a list of ids that match the search in given database
        """
        if not db:
            db = self.db
        self.ec._qs.api_key = "8fa634394aad95b85063646ab64bef2c2408"
        esr = self.ec.esearch(db=db, term=term, retmax=retmax)
        return esr.ids

    def get_paper(self, id_pubmed):
        txt = f"efetch.fcgi?db=pubmed&id={id_pubmed}&rettype=xml.xml"
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def get_ids_information(self, db_id, db=None):
        """
        This function will take  and look for the complete information of that
        for multiples Ids
        :param db: The database where we are going to look
        :param db_id: The id of the element in that database (maybe table???)
        :return:
        """
        information = []
        if not db:
            db = self.db
        print(f"{len(db_id)} papers found")
        count = 0
        total = len(db_id)
        jump = total / 100
        finished = 0
        for db_aid in db_id:
            count = count + 1

            if count >= jump:
                percent = finished/total
                print(f"Found info of {percent}%")
                count = 0

            if False:
                pass
            else:
                egs = self.ec.efetch(db=db, id=db_aid)
                cites = self.ec.elink(db=db, id=db_aid)
                if db == 'pubmed':
                    data = self.parse_pubmed_data(egs, db_aid, cites)
                    information.append(data)
                else:
                    information = getattr(egs, self.table)
            finished = finished + 1
        return information

    def get_id_information(self, db_id, db=None):
        """
        This funtion will take  and look for the complete information of that
        :param db: The database where we are going to look
        :param db_id: The id of the element in that database (maybe table???)
        :return:
        """
        if not db:
            db = self.db
        egs = self.ec.efetch(db=db, id=db_id)
        information = getattr(egs, self.table)
        return information[0]

    def parse_pubmed_data(self, egs, pid, cites):
        try:
            article_data = {'Pubmed_id': pid,
                            'Year': '0000',
                            'Title': 'NA',
                            'Key_words': [],
                            'Cite_by': '0'}
            root = egs._xml_root
            document = root.find('.//PubmedArticle/MedlineCitation/Article')
            is_article = True
            if document is None:
                document = root.find('.//PubmedBookArticle/BookDocument')
                is_article = False
            if is_article:
                year = root.find('.//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
                if year is None:
                    year_c = root.find('.//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate')
                    if year_c is not None:
                        year = year_c.text.split(" ")[0]
                    else:
                        print(f"Failed to find year for {pid}")
                        year = '0000'
                else:
                    year = year.text
            else:
                year = document.find('.//Book/PubDate/Year')
                if year is not None:
                    year = year.text
                else:
                    print(f"Failed to find year for {pid}")
                    year = '0000'

            # pubmed_info = root.find("PubmedArticle")

            # citation = pubmed_info.find("MedlineCitation")
            # article = citation.find("Article")
            # year = article.find("Journal")
            # year = year.find("JournalIssue")
            # year = year.find("PubDate")
            # year_c = year.find("Year")
            # if year_c is not None:
            #     year = year_c.text
            # else:
            #     year_c = year.find("MedlineDate")
            #     year = year_c.text.split(" ")[0]
            # title = article.find("ArticleTitle").text
            if is_article:
                title = document.find('ArticleTitle')
            else:
                title = document.find('Book/BookTitle')
            if title is None:
                print(f"Failed to find title for {pid}")
                title = 'NA'
            else:
                title = title.text
            key_words = []
            if is_article:
                citation = root.find('.//PubmedArticle/MedlineCitation/Article/')
                if citation is None:
                    print(f"Failed to find key words for {pid}")
                    citation = 0
                else:
                    key_words_xml = citation.find("KeywordList")
                    if key_words_xml is not None:
                        for c in key_words_xml:
                            key_words.append(c.text)
            else:
                #todo Add key_words for books too
                print("Key words not supported for Books yet")
            root = cites
            count = 0
            for each in root.findall('.//LinkSet/LinkSetDb/Link'):
                count = count + 1
            article_data = {'Pubmed_id': pid,
                            'Year': year,
                            'Title': title,
                            'Key_words': key_words,
                            'Cite_by': count}
        except Exception as e:
            print(e)
            print(article_data)

        return article_data
