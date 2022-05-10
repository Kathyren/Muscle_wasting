import pytest

import ncbi.main as ncbi_m
import ncbi.eutilities as eut


def test_get_papers_from_NCBI():
    ncbi_m.get_papers_from_NCBI("mirna")


def test_efetch_manula():
    """
    This test is to check the manual made efetch that is a copy with a sligthly modification to the
    one in eutilz.
    :return:
    """
    eutilcito = eut.EutilsConnection(database=eut.NCBIDatabases.Pubmed)
    r = eutilcito.efetch(db='pubmed', id='35525953')
    print(r)


def test_pruebita():
    data = '<?xml version="1.0" ?>\n<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January ' \
           '2019//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">\n<PubmedArticleSet><PubmedArticle' \
           '><MedlineCitation Status="Publisher" Owner="NLM"><PMID ' \
           'Version="1">35525953</PMID><DateRevised><Year>2022</Year><Month>05</Month><Day>07</Day></DateRevised' \
           '><Article PubModel="Electronic"><Journal><ISSN IssnType="Electronic">1471-2164</ISSN><JournalIssue ' \
           'CitedMedium="Internet"><Volume>23</Volume><Issue>1</Issue><PubDate><Year>2022</Year><Month>May</Month' \
           '><Day>07</Day></PubDate></JournalIssue><Title>BMC genomics</Title><ISOAbbreviation>BMC ' \
           'Genomics</ISOAbbreviation></Journal><ArticleTitle>Small RNA sequencing and bioinformatics analysis of ' \
           'RAW264.7-derived exosomes after Mycobacterium Bovis Bacillus Calmette-Gu\xc3\xa9rin ' \
           'infection.</ArticleTitle><Pagination><MedlinePgn>355</MedlinePgn></Pagination><ELocationID EIdType="doi" ' \
           'ValidYN="Y">10.1186/s12864-022-08590-w</ELocationID><Abstract><AbstractText Label="BACKGROUND" ' \
           'NlmCategory="BACKGROUND">The mechanisms through which Mycobacterium tuberculosis evades immune ' \
           'surveillance during tuberculosis (TB) infection remain complex. Previous studies have found that ' \
           'Mycobacteria can manipulate the miRNAs of host cells to promote their survival during host-pathogen ' \
           'interactions, and most of these effects occur at the cellular miRNA level. We attempted to investigate ' \
           'the possible related mechanisms at the exosomal miRNA level.</AbstractText><AbstractText Label="RESULTS" ' \
           'NlmCategory="RESULTS">High-throughput sequencing revealed that Bacillus Calmette-Gu\xc3\xa9rin (BCG) ' \
           'infection could alter the composition of the macrophage exosome content, and the expression levels of ' \
           'miRNAs in exosomes derived from the cell culture media of macrophages showed significant differences ' \
           'between the BCG-infected and non-infected groups. Compared with the non-infected group, 20 exosomal ' \
           'miRNAs were up-regulated and 7 exosomal miRNAs were down-regulated in the infection group (' \
           'p\xe2\x80\x89&lt;\xe2\x80\x890.05), of which mmu-miR-27b-3p, mmu-miR-93-5p, mmu-miR-25-3p, ' \
           'mmu-miR-1198-5p, mmu-let-7c-5p and let-7a-5p were significantly up-regulated. A bioinformatic analysis ' \
           'indicated that these differentially expressed exosomal miRNAs were involved in multiple biological ' \
           'processes and pathways. The target genes of top six miRNAs in up-regulated groups were positively ' \
           'correlated with the regulation of apoptosis.</AbstractText><AbstractText Label="CONCLUSIONS" ' \
           'NlmCategory="CONCLUSIONS">The expression profile of miRNA in exosomes derived from macrophage were ' \
           'altered after Mycobacterium Bovis Bacillus Calmette-Gu\xc3\xa9rin infection, and the differentially ' \
           'expressed miRNAs were involved in multiple biological processes and signalling pathways. The top six ' \
           'up-regulated miRNAs and their targeted genes were predominantly correlated with the regulation of ' \
           'apoptosis.</AbstractText><CopyrightInformation>\xc2\xa9 2022. The Author(' \
           's).</CopyrightInformation></Abstract><AuthorList CompleteYN="Y"><Author ValidYN="Y" ' \
           'EqualContrib="Y"><LastName>Zhan</LastName><ForeName>Xuehua</ForeName><Initials>X</Initials' \
           '><AffiliationInfo><Affiliation>Department of Orthopaedics, General Hospital of Ningxia Medical ' \
           'University, Yinchuan, 750004, Ningxia, China.</Affiliation></AffiliationInfo></Author><Author ValidYN="Y" ' \
           'EqualContrib="Y"><LastName>Yuan</LastName><ForeName>Wenqi</ForeName><Initials>W</Initials' \
           '><AffiliationInfo><Affiliation>Department of Orthopaedics, General Hospital of Ningxia Medical ' \
           'University, Yinchuan, 750004, Ningxia, China.</Affiliation></AffiliationInfo></Author><Author ' \
           'ValidYN="Y"><LastName>Zhou</LastName><ForeName>Yueyong</ForeName><Initials>Y</Initials><AffiliationInfo' \
           '><Affiliation>Clinical Medicine School, Ningxia Medical University, Yinchuan, 750004, Ningxia, ' \
           'China.</Affiliation></AffiliationInfo></Author><Author ' \
           'ValidYN="Y"><LastName>Ma</LastName><ForeName>Rong</ForeName><Initials>R</Initials><AffiliationInfo' \
           '><Affiliation>Department of Orthopaedics, General Hospital of Ningxia Medical University, Yinchuan, ' \
           '750004, Ningxia, China.</Affiliation></AffiliationInfo></Author><Author ' \
           'ValidYN="Y"><LastName>Ge</LastName><ForeName>Zhaohui</ForeName><Initials>Z</Initials><AffiliationInfo' \
           '><Affiliation>Department of Orthopaedics, General Hospital of Ningxia Medical University, Yinchuan, ' \
           '750004, Ningxia, China. myovid@126.com.</Affiliation></AffiliationInfo></Author></AuthorList><Language' \
           '>eng</Language><GrantList CompleteYN="Y"><Grant><GrantID>8196090244</GrantID><Agency>National Natural ' \
           'Science Foundation of China</Agency><Country/></Grant></GrantList><PublicationTypeList><PublicationType ' \
           'UI="D016428">Journal Article</PublicationType></PublicationTypeList><ArticleDate ' \
           'DateType="Electronic"><Year>2022</Year><Month>05</Month><Day>07</Day></ArticleDate></Article' \
           '><MedlineJournalInfo><Country>England</Country><MedlineTA>BMC ' \
           'Genomics</MedlineTA><NlmUniqueID>100965258</NlmUniqueID><ISSNLinking>1471-2164</ISSNLinking' \
           '></MedlineJournalInfo><CitationSubset>IM</CitationSubset><KeywordList Owner="NOTNLM"><Keyword ' \
           'MajorTopicYN="N">Apoptosis</Keyword><Keyword MajorTopicYN="N">BCG</Keyword><Keyword ' \
           'MajorTopicYN="N">Exosomal miRNA</Keyword><Keyword MajorTopicYN="N">Macrophages</Keyword><Keyword ' \
           'MajorTopicYN="N">RNA-seq</Keyword><Keyword ' \
           'MajorTopicYN="N">Tuberculosis</Keyword></KeywordList></MedlineCitation><PubmedData><History' \
           '><PubMedPubDate PubStatus="received"><Year>2022</Year><Month>01</Month><Day>26</Day></PubMedPubDate' \
           '><PubMedPubDate PubStatus="accepted"><Year>2022</Year><Month>04</Month><Day>25</Day></PubMedPubDate' \
           '><PubMedPubDate PubStatus="entrez"><Year>2022</Year><Month>5</Month><Day>7</Day><Hour>23</Hour><Minute>35' \
           '</Minute></PubMedPubDate><PubMedPubDate ' \
           'PubStatus="pubmed"><Year>2022</Year><Month>5</Month><Day>8</Day><Hour>6</Hour><Minute>0</Minute' \
           '></PubMedPubDate><PubMedPubDate PubStatus="medline"><Year>2022</Year><Month>5</Month><Day>8</Day><Hour>6' \
           '</Hour><Minute>0</Minute></PubMedPubDate></History><PublicationStatus>epublish</PublicationStatus' \
           '><ArticleIdList><ArticleId IdType="pubmed">35525953</ArticleId><ArticleId ' \
           'IdType="doi">10.1186/s12864-022-08590-w</ArticleId><ArticleId ' \
           'IdType="pii">10.1186/s12864-022-08590-w</ArticleId></ArticleIdList><ReferenceList><Reference><Citation' \
           '>World Health Organization. Global Tuberculosis Report 2021. Geneva: World Health Organization; 2021. p. ' \
           '2021.</Citation></Reference><Reference><Citation>Liu CH, Liu H, Ge B. Innate immunity in tuberculosis: ' \
           'host defense vs pathogen evasion. Cell Mol Immunol. 2017;14(' \
           '12):963\xe2\x80\x9375.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1038/cmi.2017.88</ArticleId></ArticleIdList></Reference><Reference><Citation>Zhai W, ' \
           'Wu F, Zhang Y, Fu Y, Liu Z. The immune escape mechanisms of mycobacterium tuberculosis. Int J Mol Sci. ' \
           '2019;20(2):340.</Citation></Reference><Reference><Citation>Abdalla AE, Ejaz H, Mahjoob MO, Alameen AAM, ' \
           'Abosalif KOA, Elamir MYM, et al. Intelligent mechanisms of macrophage apoptosis subversion by ' \
           'mycobacterium. Pathogens. 2020;9(3):218.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.3390/pathogens9030218</ArticleId></ArticleIdList></Reference><Reference><Citation>Alipoor ' \
           'SD, Tabarsi P, Varahram M, Movassaghi M, Dizaji MK, Folkerts G, et al. Serum exosomal miRNAs are ' \
           'associated with active pulmonary tuberculosis. Dis Markers. ' \
           '2019;2019:1907426.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1155/2019/1907426</ArticleId></ArticleIdList></Reference><Reference><Citation>Ha M, ' \
           'Kim VN. Regulation of microRNA biogenesis. Nat Rev Mol Cell Biol. 2014;15(' \
           '8):509\xe2\x80\x9324.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1038/nrm3838</ArticleId></ArticleIdList></Reference><Reference><Citation>Jadli AS, ' \
           'Ballasy N, Edalat P, Patel VB. Inside (sight) of tiny communicator: exosome biogenesis, secretion, ' \
           'and uptake. Mol Cell Biochem. 2020;467(' \
           '1\xe2\x80\x932):77\xe2\x80\x9394.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1007/s11010-020-03703-z</ArticleId></ArticleIdList></Reference><Reference><Citation>O' \
           '\'Brien K, Breyne K, Ughetto S, Laurent LC, Breakefield XO. RNA delivery by extracellular vesicles in ' \
           'mammalian cells and its applications. Nat Rev MolCell Biol. 2020;21(' \
           '10):585\xe2\x80\x93606.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1038/s41580-020-0251-y</ArticleId></ArticleIdList></Reference><Reference><Citation>Chen ' \
           'WX, Liu XM, Lv MM, Chen L, Zhao JH, Zhong SL, et al. Exosomes from drug-resistant breast cancer cells ' \
           'transmit chemoresistance by a horizontal transfer of microRNAs. PLoS One. 2014;9(' \
           '4):e95240.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1371/journal.pone.0095240</ArticleId></ArticleIdList></Reference><Reference><Citation' \
           '>Ismail N, Wang Y, Dakhlallah D, Moldovan L, Agarwal K, Batte K, et al. Macrophage microvesicles induce ' \
           'macrophage differentiation and miR-223 transfer. Blood. 2013;121(' \
           '6):984\xe2\x80\x9395.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1182/blood-2011-08-374793</ArticleId></ArticleIdList></Reference><Reference><Citation' \
           '>Singh PP, Li L, Schorey JS. Exosomal RNA from Mycobacterium tuberculosis-infected cells is functional in ' \
           'recipient macrophages. Traffic. 2015;16(6):555\xe2\x80\x9371.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1111/tra.12278</ArticleId></ArticleIdList></Reference><Reference><Citation>Bettencourt P, ' \
           'Carmo N, Pires D, Tim\xc3\xb3teo P, Anes E. Mycobacterial infection of macrophages: the effect of the ' \
           'multiplicity of infection. In: M\xc3\xa9ndez-Vilas A, editor. Antimicrobial research: novel bioknowledge ' \
           'and educational programs. Badajoz: Formatex Research Center; 2017. p. ' \
           '651\xe2\x80\x9364.</Citation></Reference><Reference><Citation>Szklarczyk D, Gable AL, Lyon D, Junge A, ' \
           'Wyder S, Huerta-Cepas J, et al. STRING v11: protein-protein association networks with increased coverage, ' \
           'supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Res. 2019;47(' \
           'D1):D607\xe2\x80\x93D13.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1093/nar/gky1131</ArticleId></ArticleIdList></Reference><Reference><Citation>Stark C, ' \
           'Breitkreutz BJ, Reguly T, Boucher L, Breitkreutz A, Tyers M. BioGRID: a general repository for ' \
           'interaction datasets. Nucleic Acids Res. 2006;34(Database ' \
           'issue):D535\xe2\x80\x939.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1093/nar/gkj109</ArticleId></ArticleIdList></Reference><Reference><Citation>Li T, ' \
           'Wernersson R, Hansen RB, Horn H, Mercer J, Slodkowicz G, et al. A scored human protein-protein ' \
           'interaction network to catalyze genomic interpretation. Nat Methods. 2017;14(' \
           '1):61\xe2\x80\x934.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1038/nmeth.4083</ArticleId></ArticleIdList></Reference><Reference><Citation>Bader GD, ' \
           'Hogue CW. An automated method for finding molecular complexes in large protein interaction networks. BMC ' \
           'Bioinform. 2003;4:2.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1186/1471-2105-4-2</ArticleId></ArticleIdList></Reference><Reference><Citation>Shannon P, ' \
           'Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, et al. Cytoscape: a software environment for integrated ' \
           'models of biomolecular interaction networks. Genome Res. 2003;13(' \
           '11):2498\xe2\x80\x93504.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1101/gr.1239303</ArticleId></ArticleIdList></Reference><Reference><Citation>Singh PP, ' \
           'Smith VL, Karakousis PC, Schorey JS. Exosomes isolated from mycobacteria-infected mice or cultured ' \
           'macrophages can recruit and activate immune cells in vitro and in vivo. J Immunol. 2012;189(' \
           '2):777\xe2\x80\x9385.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.4049/jimmunol.1103638</ArticleId></ArticleIdList></Reference><Reference><Citation>Alipoor ' \
           'SD, Adcock IM, Garssen J, Mortaz E, Varahram M, Mirsaeidi M, et al. The roles of miRNAs as potential ' \
           'biomarkers in lung diseases. Eur J Pharmacol. ' \
           '2016;791:395\xe2\x80\x93404.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1016/j.ejphar.2016.09.015</ArticleId></ArticleIdList></Reference><Reference><Citation' \
           '>Behar SM, Martin CJ, Booty MG, Nishimura T, Zhao X, Gan HX, et al. Apoptosis is an innate defense ' \
           'function of macrophages against mycobacterium tuberculosis. Mucosal Immunol. 2011;4(' \
           '3):279\xe2\x80\x9387.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1038/mi.2011.3</ArticleId></ArticleIdList></Reference><Reference><Citation>Alipoor SD, ' \
           'Adcock IM, Folkerts G, Garssen J, Mortaz E. A bioinformatics analysis of exosomal microRNAs released ' \
           'following mycobacterial infection. Int J Mycobacteriol. 2019;8(' \
           '3):218\xe2\x80\x9322.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.4103/ijmy.ijmy_88_19</ArticleId></ArticleIdList></Reference><Reference><Citation>Agarwal ' \
           'RG, Sharma P, Nyati KK. microRNAs in mycobacterial infection: modulation of host immune response and ' \
           'apoptotic pathways. Immune Netw. 2019;19(5):e30.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.4110/in.2019.19.e30</ArticleId></ArticleIdList></Reference><Reference><Citation>Kim JK, ' \
           'Kim TS, Basu J, Jo EK. MicroRNA in innate immunity and autophagy during mycobacterial infection. Cell ' \
           'Microbiol. 2017;19(1):e12687.</Citation></Reference><Reference><Citation>Loeuillet C, Martinon F, ' \
           'Perez C, Munoz M, Thome M, Meylan PR. Mycobacterium tuberculosis subverts innate immunity to evade ' \
           'specific effectors. J Immunol. 2006;177(9):6245\xe2\x80\x9355.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.4049/jimmunol.177.9.6245</ArticleId></ArticleIdList></Reference><Reference><Citation>Fu ' \
           'Y, Yi Z, Wu X, Li J, Xu F. Circulating microRNAs in patients with active pulmonary tuberculosis. J Clin ' \
           'Microbiol. 2011;49(12):4246\xe2\x80\x9351.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1128/JCM.05459-11</ArticleId></ArticleIdList></Reference><Reference><Citation>Sharbati J, ' \
           'Lewin A, Kutz-Lohroff B, Kamal E, Einspanier R, Sharbati S. Integrated microRNA-mRNA-analysis of human ' \
           'monocyte derived macrophages upon mycobacterium avium subsp. hominissuis infection. PLoS One. 2011;6(' \
           '5):e20258.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1371/journal.pone.0020258</ArticleId></ArticleIdList></Reference><Reference><Citation' \
           '>Bettencourt P, Marion S, Pires D, Santos LF, Lastrucci C, Carmo N, et al. Actin-binding protein ' \
           'regulation by microRNAs as a novel microbial strategy to modulate phagocytosis by host cells: the case of ' \
           'N-wasp and miR-142-3p. Front Cell Infect Microbiol. 2013;3:19.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.3389/fcimb.2013.00019</ArticleId></ArticleIdList></Reference><Reference><Citation>Guo L, ' \
           'Zhou L, Gao Q, Zhang A, Wei J, Hong D, et al. MicroRNA-144-3p inhibits autophagy activation and enhances ' \
           'Bacillus Calmette-Guerin infection by targeting ATG4a in RAW264.7 macrophage cells. PLoS One. 2017;12(' \
           '6):e0179772.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1371/journal.pone.0179772</ArticleId></ArticleIdList></Reference><Reference><Citation>Gu ' \
           'X, Gao Y, Mu DG, Fu EQ. MiR-23a-5p modulates mycobacterial survival and autophagy during mycobacterium ' \
           'tuberculosis infection through TLR2/MyD88/NF-kappaB pathway by targeting TLR2. Exp Cell Res. 2017;354(' \
           '2):71\xe2\x80\x937.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1016/j.yexcr.2017.03.039</ArticleId></ArticleIdList></Reference><Reference><Citation' \
           '>Alipoor SD, Mortaz E, Tabarsi P, Farnia P, Mirsaeidi M, Garssen J, et al. Bovis bacillus calmette-guerin ' \
           '(BCG) infection induces exosomal miRNA release by human macrophages. J Transl Med. 2017;15(' \
           '1):105.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1186/s12967-017-1205-9</ArticleId></ArticleIdList></Reference><Reference><Citation' \
           '>Alipoor SD, Mortaz E, Tabarsi P, Marjani M, Varahram M, Folkerts G, et al. miR-1224 expression is ' \
           'increased in human macrophages after infection with bacillus calmette-guerin (BCG). Iran J Allergy Asthma ' \
           'Immunol. 2018;17(3):250\xe2\x80\x937.</Citation><ArticleIdList><ArticleId ' \
           'IdType="pubmed">29908542</ArticleId></ArticleIdList></Reference><Reference><Citation>Daniel J, Maamar H, ' \
           'Deb C, Sirakova TD, Kolattukudy PE. Mycobacterium tuberculosis uses host triacylglycerol to accumulate ' \
           'lipid droplets and acquires a dormancy-like phenotype in lipid-loaded macrophages. PLoS Pathog. 2011;7(' \
           '6):e1002093.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1371/journal.ppat.1002093</ArticleId></ArticleIdList></Reference><Reference><Citation>Kim ' \
           'MJ, Wainwright HC, Locketz M, Bekker LG, Walther GB, Dittrich C, et al. Caseation of human tuberculosis ' \
           'granulomas correlates with elevated host lipid metabolism. EMBO Mol Med. 2010;2(' \
           '7):258\xe2\x80\x9374.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1002/emmm.201000079</ArticleId></ArticleIdList></Reference><Reference><Citation>Zheng X, ' \
           'Ye C, Zhao J, Bian P, Zhang Y, Jia Z. Alterations and clinical signifecance of exosome-containing innate ' \
           'immunity related lncRNAs in patients of hemorrhagic fever with renal syndrome. Xi Bao Yu Fen Zi Mian Yi ' \
           'Xue Za Zhi. 2016;32(11):1522\xe2\x80\x936.</Citation><ArticleIdList><ArticleId ' \
           'IdType="pubmed">27774948</ArticleId></ArticleIdList></Reference><Reference><Citation>Zhang D, Yi Z, ' \
           'Fu Y. Downregulation of miR-20b-5p facilitates Mycobacterium tuberculosis survival in RAW 264.7 ' \
           'macrophages via attenuating the cell apoptosis by Mcl-1 upregulation. J Cell Biochem. 2019;120(' \
           '4):5889\xe2\x80\x9396.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1002/jcb.27874</ArticleId></ArticleIdList></Reference><Reference><Citation>Hadifar S, ' \
           'Fateh A, Yousefi MH, Siadat SD, Vaziri F. Exosomes in tuberculosis: Still terra incognita? J Cell ' \
           'Physiol. 2019;234(3):2104\xe2\x80\x9311.</Citation><ArticleIdList><ArticleId ' \
           'IdType="doi">10.1002/jcp.27555</ArticleId></ArticleIdList></Reference></ReferenceList></PubmedData' \
           '></PubmedArticle></PubmedArticleSet> '
    import xml.etree.ElementTree as ET
    root = ET.fromstring(data)
    pubmed_info = root.find("PubmedArticle")
    citation = pubmed_info.find("MedlineCitation")
    article = citation.find("Article")
    year = article.find("Journal")
    year = year.find("JournalIssue")
    year = year.find("PubDate")
    year = year.find("Year").text
    title = article.find("ArticleTitle").text
    key_words_xml = citation.find("KeywordList")
    key_words = []
    if key_words_xml is not None:
        for c in key_words_xml:
            key_words.append(c.text)
            print(c)
    x = root.find("PubmedData")
    x = x


def test_get_ids_information():
    eutilcito = eut.EutilsConnection(database=eut.NCBIDatabases.Pubmed)
    eutilcito.get_ids_information(db_id=[ 35521437, 35524416, 21876726], db='pubmed')
