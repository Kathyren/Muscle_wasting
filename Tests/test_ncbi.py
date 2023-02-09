import pytest
from eutils import Client

import ncbi.main as ncbi_m
import ncbi.eutilities as eut
import xml.etree.ElementTree as ET

eutilcito = eut.EutilsConnection(database=eut.NCBIDatabases.Pubmed)

book_info_one_xml = b'<?xml version="1.0" ?>\n<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January ' \
                    b'2019//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">\n<PubmedArticleSet' \
                    b'><PubmedBookArticle><BookDocument><PMID Version="1">34662018</PMID><ArticleIdList><ArticleId ' \
                    b'IdType="bookaccession">NBK574504</ArticleId></ArticleIdList><Book><Publisher><PublisherName' \
                    b'>StatPearls Publishing</PublisherName><PublisherLocation>Treasure Island (' \
                    b'FL)</PublisherLocation></Publisher><BookTitle ' \
                    b'book="statpearls">StatPearls</BookTitle><PubDate><Year>2022</Year><Month>01</Month></PubDate' \
                    b'><BeginningDate><Year>2022</Year><Month>01</Month></BeginningDate><Medium>Internet</Medium></Book' \
                    b'><ArticleTitle book="statpearls" part="article-139270">Fertility Sparing Management In Uterine ' \
                    b'Fibroids</ArticleTitle><Language>eng</Language><AuthorList ' \
                    b'Type="authors"><Author><LastName>Rezk</LastName><ForeName>Andrew</ForeName><Initials>A</Initials' \
                    b'><AffiliationInfo><Affiliation>University of Miami Leonard M. Miller School of Medicine,' \
                    b'</Affiliation></AffiliationInfo></Author><Author><LastName>Kahn</LastName><ForeName>Jenna' \
                    b'</ForeName><Initials>J</Initials><AffiliationInfo><Affiliation>Montefiore Medical Center/ Albert ' \
                    b'Einstein College of Medicine</Affiliation></AffiliationInfo></Author><Author><LastName>Singh' \
                    b'</LastName><ForeName>Manvinder</ForeName><Initials>M</Initials><AffiliationInfo><Affiliation' \
                    b'>Montefiore/Einstein</Affiliation></AffiliationInfo></Author></AuthorList><PublicationType ' \
                    b'UI="D000072643">Study Guide</PublicationType><Abstract><AbstractText>Uterine leiomyoma (fibroids) ' \
                    b'are the most common benign gynecologic tumors, occurring in up to 70% of women by menopause.[' \
                    b'1]\xc2\xa0While many fibroids are asymptomatic and are only discovered incidentally, 25 to 30 % of ' \
                    b'women experience a spectrum of symptoms that increase morbidity and adversely affect their quality ' \
                    b'of life.[2]\xc2\xa0 The most common symptoms include abnormal uterine bleeding (AUB), ' \
                    b'heavy menstrual bleeding (HMB), pelvic pain and pressure, anemia, and bladder and/or bowel ' \
                    b'dysfunction. Importantly, the presence of one or multiple fibroids may affect fertility, ' \
                    b'as the distortion of the uterus can prevent successful implantation and/or continued survival of ' \
                    b'an intrauterine pregnancy.\xc2\xa0 The presence of fibroids can introduce obstetrical ' \
                    b'complications, such as recurrent pregnancy loss (RPL), preterm labor (PTL), abnormal placentation, ' \
                    b'increased rates of cesarean section, and postpartum hemorrhage.\xc2\xa0 With the current trend of ' \
                    b'increased childbearing age, fertility-sparing management of uterine fibroids is critical.\xc2\xa0 ' \
                    b'In this review, we discuss the clinical presentation of fibroids in women of reproductive age and ' \
                    b'the medical and surgical management options available for women with fibroids who wish to become ' \
                    b'pregnant in the future.\xc2\xa0</AbstractText><CopyrightInformation>Copyright \xc2\xa9 2022, ' \
                    b'StatPearls Publishing LLC.</CopyrightInformation></Abstract><Sections><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" sec="article-139270.s1">Continuing Education ' \
                    b'Activity</SectionTitle></Section><Section><SectionTitle book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s2">Introduction</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s3">Etiology</SectionTitle></Section><Section><SectionTitle book="statpearls" ' \
                    b'part="article-139270" sec="article-139270.s4">Epidemiology</SectionTitle></Section><Section' \
                    b'><SectionTitle book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s5">Pathophysiology</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s6">Histopathology</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" sec="article-139270.s7">History and ' \
                    b'Physical</SectionTitle></Section><Section><SectionTitle book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s8">Evaluation</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" sec="article-139270.s9">Treatment / ' \
                    b'Management</SectionTitle></Section><Section><SectionTitle book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s10">Differential Diagnosis</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s11">Prognosis</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s12">Complications</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" sec="article-139270.s13">Deterrence and Patient ' \
                    b'Education</SectionTitle></Section><Section><SectionTitle book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s14">Enhancing Healthcare Team Outcomes ' \
                    b'</SectionTitle></Section><Section><SectionTitle book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s15">Review Questions</SectionTitle></Section><Section><SectionTitle ' \
                    b'book="statpearls" part="article-139270" ' \
                    b'sec="article-139270.s21">References</SectionTitle></Section></Sections><ContributionDate><Year' \
                    b'>2021</Year><Month>11</Month><Day>4</Day></ContributionDate></BookDocument><PubmedBookData' \
                    b'><History><PubMedPubDate PubStatus="pubmed"><Year>2021</Year><Month>10</Month><Day>19</Day><Hour>6' \
                    b'</Hour><Minute>1</Minute></PubMedPubDate><PubMedPubDate ' \
                    b'PubStatus="medline"><Year>2021</Year><Month>10</Month><Day>19</Day><Hour>6</Hour><Minute>1</Minute' \
                    b'></PubMedPubDate><PubMedPubDate ' \
                    b'PubStatus="entrez"><Year>2021</Year><Month>10</Month><Day>19</Day><Hour>6</Hour><Minute>1</Minute' \
                    b'></PubMedPubDate></History><PublicationStatus>ppublish</PublicationStatus><ArticleIdList' \
                    b'><ArticleId IdType="pubmed">34662018</ArticleId></ArticleIdList></PubmedBookData' \
                    b'></PubmedBookArticle></PubmedArticleSet> '
no_cites_one_xml = b'<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD elink 20101123//EN" ' \
                   b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n<eLinkResult>\n\n  <LinkSet>\n    ' \
                   b'<DbFrom>pubmed</DbFrom>\n    <IdList>\n      <Id>34662018</Id>\n    </IdList>\n    \n  ' \
                   b'</LinkSet>\n</eLinkResult>\n '
tons_cites_one_xml = b'<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD elink ' \
                     b'20101123//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n<eLinkResult>\n\n  ' \
                     b'<LinkSet>\n    <DbFrom>pubmed</DbFrom>\n    <IdList>\n      <Id>21876726</Id>\n    </IdList>\n    ' \
                     b'<LinkSetDb>\n      <DbTo>pmc</DbTo>\n      <LinkName>pubmed_pmc_refs</LinkName>\n      \n        ' \
                     b'<Link>\n\t\t\t\t<Id>9068894</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>9013908</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8991685</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8980040</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8906861</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8888139</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8876091</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8875259</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8865617</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8835344</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8819181</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8810016</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8770970</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8766409</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8722215</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8721122</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8706681</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8667074</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8657631</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8637759</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8584775</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8559595</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8537887</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8526851</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8500979</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8455141</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8440450</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8431789</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8417445</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8407884</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8407881</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8345803</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8327470</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8323226</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8320674</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8297715</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8286722</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8277821</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8176998</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8165805</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8160712</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8131927</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8104913</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8071896</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8025741</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>8023533</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7998656</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7994852</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7958126</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7944561</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7944388</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7906989</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7905690</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7901554</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7887693</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7840842</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7814576</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7783313</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7774862</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7753098</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7744002</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7734763</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7695640</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7692229</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7595554</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7578871</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7564998</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7556308</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7523352</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7519665</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7496643</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7431678</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7414114</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7412956</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7409063</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7400702</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7391522</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7378182</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7362528</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7347774</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7345075</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7340524</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7318537</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7313400</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7304098</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7302792</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7283084</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7277307</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7227815</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7226815</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7226573</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7225259</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7217362</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7205224</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7196212</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7186296</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7183550</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7181663</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7145962</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7099601</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7090279</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7081927</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7072607</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7068154</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7048814</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>7016778</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6996099</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6991626</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6942270</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6918866</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6912723</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6902397</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6900535</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6888613</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6872767</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6858646</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6854698</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6843021</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6831938</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6828747</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6828701</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6815838</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6794619</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6785768</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6784011</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6780992</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6753017</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6745844</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6739429</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6708438</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6694442</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6688994</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6627924</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6627787</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6609620</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6604862</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6591556</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6566274</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6548868</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6541292</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6491393</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6488179</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6454737</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6433537</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6428090</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6394792</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6386845</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6381937</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6359271</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6349164</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6345809</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6332863</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6325302</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6317223</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6305133</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6302004</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6299201</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6277517</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6251105</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6249419</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6234778</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6205447</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6205207</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6198229</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6187740</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6182003</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6164282</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6150479</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6096814</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6092208</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6089553</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6083546</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6055191</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6048473</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6045325</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6038041</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6036139</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6032513</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6029891</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>6008571</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5991072</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5940186</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5922373</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5912192</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5903829</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5887919</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5880298</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5859015</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5855813</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5842184</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5836374</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5835801</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5818524</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5809724</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5804051</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5801828</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5790938</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5783023</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5780426</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5754134</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5745058</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5717444</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5712011</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5706431</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5688943</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5666769</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5662033</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5660123</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5648287</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5640880</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5621885</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5618269</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5610870</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5599553</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5579732</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5570485</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5548323</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5548308</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5537487</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5522217</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5519199</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5517456</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5489338</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5466984</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5463428</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5454934</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5426696</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5420547</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5390430</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5370134</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5363297</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5353755</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5341331</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5309960</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5291278</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5244794</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5222964</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5154766</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5147504</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5131095</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5114451</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5110222</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5101684</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5066681</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5047747</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5047620</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>5035490</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4994100</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4992611</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4991373</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4990609</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4980326</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4967388</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4966537</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4964959</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4949450</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4947197</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4921266</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4913210</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4911875</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4888965</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4888503</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4887680</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4881453</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4875887</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4867787</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4859755</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4859494</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4853860</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4840335</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4834364</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4811151</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4797991</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4773778</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4773204</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4765724</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4745815</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4716143</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4709998</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4692419</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4685040</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4675626</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4656109</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4655533</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4642616</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4637461</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4625016</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4594673</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4587316</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4578249</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4578234</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4577581</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4572411</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4550760</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4514545</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4500606</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4496044</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4463923</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4454760</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4449452</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4448148</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4443736</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4435734</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4423606</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4412282</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4379073</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4370721</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4333915</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4330530</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4327859</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4321981</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4317028</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4261509</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4261240</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4258658</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4254633</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4180193</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4176915</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4154666</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4150742</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4148874</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4145399</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4128075</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4104030</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4102923</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4082760</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4073282</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4063995</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4063195</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4063074</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4059258</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4057926</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4047354</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4044516</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4035236</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4030698</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>4005017</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3982321</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3943742</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3930942</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3928468</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3917938</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3899882</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3896203</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3883893</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3879572</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3873329</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3848851</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3835174</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3830597</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3823269</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3817212</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3812991</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3804457</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3782363</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3624465</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3591653</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3580451</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3580056</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3560214</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3511851</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3498077</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3496552</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3495063</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3469598</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3433256</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3422911</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3355967</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3351299</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3340916</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3339323</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3330238</Id>\n\t\t\t</Link>\n        ' \
                     b'<Link>\n\t\t\t\t<Id>3235036</Id>\n\t\t\t</Link>\n      \n    </LinkSetDb>\n  ' \
                     b'</LinkSet>\n</eLinkResult>\n '

cites_many_xml = b'<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD elink 20101123//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n<eLinkResult>\n\n  <LinkSet>\n    <DbFrom>pubmed</DbFrom>\n    <IdList>\n      <Id>21876726</Id>\n    </IdList>\n    <LinkSetDb>\n      <DbTo>pmc</DbTo>\n      <LinkName>pubmed_pmc_refs</LinkName>\n      \n        <Link>\n\t\t\t\t<Id>9094574</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>9068894</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>9013908</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8991685</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8980040</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8906861</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8888139</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8876091</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8875259</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8865617</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8835344</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8819181</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8810016</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8770970</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8766409</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8722215</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8721122</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8706681</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8667074</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8657631</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8637759</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8584775</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8559595</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8537887</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8526851</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8500979</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8455141</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8440450</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8431789</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8417445</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8407884</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8407881</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8345803</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8327470</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8323226</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8320674</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8297715</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8286722</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8277821</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8176998</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8165805</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8160712</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8131927</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8104913</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8071896</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8025741</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8023533</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7998656</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7994852</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7958126</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7944561</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7944388</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7906989</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7905690</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7901554</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7887693</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7840842</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7814576</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7783313</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7774862</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7753098</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7744002</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7734763</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7695640</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7692229</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7595554</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7578871</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7564998</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7556308</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7523352</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7519665</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7496643</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7431678</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7414114</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7412956</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7409063</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7400702</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7391522</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7378182</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7362528</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7347774</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7345075</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7340524</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7318537</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7313400</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7304098</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7302792</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7283084</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7277307</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7227815</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7226815</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7226573</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7225259</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7217362</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7205224</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7196212</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7186296</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7183550</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7181663</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7145962</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7099601</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7090279</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7081927</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7072607</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7068154</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7048814</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7016778</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6996099</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6991626</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6942270</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6918866</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6912723</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6902397</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6900535</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6888613</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6872767</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6858646</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6854698</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6843021</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6831938</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6828747</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6828701</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6815838</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6794619</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6785768</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6784011</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6780992</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6753017</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6745844</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6739429</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6708438</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6694442</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6688994</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6627924</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6627787</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6609620</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6604862</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6591556</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6566274</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6548868</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6541292</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6491393</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6488179</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6454737</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6433537</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6428090</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6394792</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6386845</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6381937</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6359271</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6349164</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6345809</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6332863</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6325302</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6317223</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6305133</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6302004</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6299201</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6277517</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6251105</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6249419</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6234778</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6205447</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6205207</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6198229</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6187740</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6182003</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6164282</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6150479</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6096814</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6092208</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6089553</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6083546</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6055191</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6048473</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6045325</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6038041</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6036139</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6032513</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6029891</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6008571</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5991072</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5940186</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5922373</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5912192</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5903829</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5887919</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5880298</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5859015</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5855813</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5842184</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5836374</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5835801</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5818524</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5809724</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5804051</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5801828</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5790938</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5783023</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5780426</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5754134</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5745058</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5717444</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5712011</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5706431</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5688943</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5666769</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5662033</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5660123</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5648287</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5640880</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5621885</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5618269</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5610870</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5599553</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5579732</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5570485</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5548323</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5548308</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5537487</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5522217</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5519199</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5517456</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5489338</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5466984</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5463428</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5454934</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5426696</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5420547</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5390430</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5370134</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5363297</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5353755</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5341331</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5309960</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5291278</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5244794</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5222964</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5154766</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5147504</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5131095</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5114451</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5110222</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5101684</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5066681</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5047747</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5047620</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5035490</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4994100</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4992611</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4991373</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4990609</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4980326</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4967388</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4966537</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4964959</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4949450</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4947197</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4921266</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4913210</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4911875</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4888965</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4888503</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4887680</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4881453</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4875887</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4867787</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4859755</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4859494</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4853860</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4840335</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4834364</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4811151</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4797991</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4773778</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4773204</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4765724</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4745815</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4716143</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4709998</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4692419</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4685040</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4675626</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4656109</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4655533</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4642616</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4637461</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4625016</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4594673</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4587316</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4578249</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4578234</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4577581</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4572411</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4550760</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4514545</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4500606</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4496044</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4463923</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4454760</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4449452</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4448148</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4443736</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4435734</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4423606</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4412282</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4379073</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4370721</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4333915</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4330530</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4327859</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4321981</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4317028</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4261509</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4261240</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4258658</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4254633</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4180193</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4176915</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4154666</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4150742</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4148874</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4145399</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4128075</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4104030</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4102923</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4082760</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4073282</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4063995</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4063195</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4063074</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4059258</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4057926</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4047354</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4044516</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4035236</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4030698</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4005017</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3982321</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3943742</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3930942</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3928468</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3917938</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3899882</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3896203</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3883893</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3879572</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3873329</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3848851</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3835174</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3830597</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3823269</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3817212</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3812991</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3804457</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3782363</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3624465</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3591653</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3580451</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3580056</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3560214</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3511851</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3498077</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3496552</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3495063</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3469598</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3433256</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3422911</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3355967</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3351299</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3340916</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3339323</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3330238</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3235036</Id>\n\t\t\t</Link>\n      \n    </LinkSetDb>\n  </LinkSet>\n\n  <LinkSet>\n    <DbFrom>pubmed</DbFrom>\n    <IdList>\n      <Id>21876761</Id>\n    </IdList>\n    <LinkSetDb>\n      <DbTo>pmc</DbTo>\n      <LinkName>pubmed_pmc_refs</LinkName>\n      \n        <Link>\n\t\t\t\t<Id>9030904</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8874453</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8749911</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8369991</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8308526</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>8092845</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7922494</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7904927</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7793277</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7721145</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>7332439</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6893788</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>6556556</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5923461</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5877976</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5800646</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5712301</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5642910</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5577925</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5545684</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5364581</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5331900</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5302377</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5133561</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5124249</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>5011806</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4862634</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4793885</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4645647</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>4248897</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3999893</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3857596</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3809390</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3700278</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3694120</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3637789</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3624376</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3564832</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3553009</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3543106</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3446598</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3419081</Id>\n\t\t\t</Link>\n        <Link>\n\t\t\t\t<Id>3364951</Id>\n\t\t\t</Link>\n      \n    </LinkSetDb>\n  </LinkSet>\n</eLinkResult>\n'
article_info_many_xml = b'<?xml version="1.0" ?>\n<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, ' \
                        b'1st January 2019//EN" ' \
                        b'"https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">\n<PubmedArticleSet' \
                        b'><PubmedArticle><MedlineCitation Status="PubMed-not-MEDLINE" Owner="NLM"><PMID ' \
                        b'Version="1">21876726</PMID><DateCompleted><Year>2011</Year><Month>11</Month><Day>10</Day' \
                        b'></DateCompleted><DateRevised><Year>2022</Year><Month>03</Month><Day>31</Day></DateRevised' \
                        b'><Article PubModel="Print-Electronic"><Journal><ISSN ' \
                        b'IssnType="Electronic">1757-594X</ISSN><JournalIssue ' \
                        b'CitedMedium="Internet"><Volume>3</Volume><PubDate><Year>2011</Year></PubDate></JournalIssue' \
                        b'><Title>F1000 biology reports</Title><ISOAbbreviation>F1000 Biol ' \
                        b'Rep</ISOAbbreviation></Journal><ArticleTitle>Exosomes: secreted vesicles and intercellular ' \
                        b'communications.</ArticleTitle><Pagination><MedlinePgn>15</MedlinePgn></Pagination' \
                        b'><ELocationID EIdType="doi" ' \
                        b'ValidYN="Y">10.3410/B3-15</ELocationID><Abstract><AbstractText>Exosomes are small membrane ' \
                        b'vesicles of endocytic origin secreted by most cell types, and are thought to play important ' \
                        b'roles in intercellular communications. Although exosomes were originally described in 1983, ' \
                        b'interest in these vesicles has really increased dramatically in the last 3 years, ' \
                        b'after the finding that they contain mRNA and microRNA. This discovery sparked renewed ' \
                        b'interest for the general field of membrane vesicles involved in intercellular ' \
                        b'communications, and research on these structures has grown exponentially over the last few ' \
                        b'years, probing their composition and function, as well as their potential value as ' \
                        b'biomarkers.</AbstractText></Abstract><AuthorList CompleteYN="Y"><Author ' \
                        b'ValidYN="Y"><LastName>Th\xc3\xa9ry</LastName><ForeName>Clotilde</ForeName><Initials>C' \
                        b'</Initials><AffiliationInfo><Affiliation>Institut Curie INSERM U932, ' \
                        b'Paris France.</Affiliation></AffiliationInfo></Author></AuthorList><Language>eng</Language' \
                        b'><PublicationTypeList><PublicationType UI="D016428">Journal ' \
                        b'Article</PublicationType></PublicationTypeList><ArticleDate ' \
                        b'DateType="Electronic"><Year>2011</Year><Month>07</Month><Day>01</Day></ArticleDate' \
                        b'></Article><MedlineJournalInfo><Country>England</Country><MedlineTA>F1000 Biol ' \
                        b'Rep</MedlineTA><NlmUniqueID>101506835</NlmUniqueID><ISSNLinking>1757-594X</ISSNLinking' \
                        b'></MedlineJournalInfo></MedlineCitation><PubmedData><History><PubMedPubDate ' \
                        b'PubStatus="entrez"><Year>2011</Year><Month>8</Month><Day>31</Day><Hour>6</Hour><Minute>0' \
                        b'</Minute></PubMedPubDate><PubMedPubDate ' \
                        b'PubStatus="pubmed"><Year>2011</Year><Month>8</Month><Day>31</Day><Hour>6</Hour><Minute>0' \
                        b'</Minute></PubMedPubDate><PubMedPubDate ' \
                        b'PubStatus="medline"><Year>2011</Year><Month>8</Month><Day>31</Day><Hour>6</Hour><Minute>1' \
                        b'</Minute></PubMedPubDate></History><PublicationStatus>ppublish</PublicationStatus' \
                        b'><ArticleIdList><ArticleId IdType="pubmed">21876726</ArticleId><ArticleId ' \
                        b'IdType="doi">10.3410/B3-15</ArticleId><ArticleId IdType="pii">15</ArticleId><ArticleId ' \
                        b'IdType="pmc">PMC3155154</ArticleId></ArticleIdList><ReferenceList><Reference><Citation>Nat ' \
                        b'Cell Biol. 2008 Dec;10(12):1470-6</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19011622</ArticleId></ArticleIdList></Reference><Reference><Citation>Nat ' \
                        b'Cell Biol. 2010 Jan;12(1):19-30; sup pp 1-13</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19966785</ArticleId></ArticleIdList></Reference><Reference><Citation>Nat ' \
                        b'Med. 1998 May;4(5):594-600</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">9585234</ArticleId></ArticleIdList></Reference><Reference><Citation>J Cell ' \
                        b'Biol. 2010 Apr 19;189(2):223-32</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">20404108</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Immunol. 1999 Oct 15;163(8):4564-73</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">10510400</ArticleId></ArticleIdList></Reference><Reference><Citation>J Cell ' \
                        b'Biol. 1983 Aug;97(2):329-39</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">6309857</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Immunol. 2005 Aug 15;175(4):2237-43</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16081791</ArticleId></ArticleIdList></Reference><Reference><Citation>Stem ' \
                        b'Cell Rev. 2008 Sep;4(3):137-47</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">18665337</ArticleId></ArticleIdList></Reference><Reference><Citation>Nat ' \
                        b'Immunol. 2002 Dec;3(12):1156-62</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">12426563</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 2006 Jul 25;103(30):11172-7</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16837572</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Traffic. 2005 Feb;6(2):131-43</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">15634213</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 2004 Jun 29;101(26):9683-8</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">15210972</ArticleId></ArticleIdList></Reference><Reference><Citation>Nat ' \
                        b'Med. 2001 Mar;7(3):297-303</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">11231627</ArticleId></ArticleIdList></Reference><Reference><Citation>Stem ' \
                        b'Cell Res. 2010 May;4(3):214-22</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">20138817</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Science. 2008 Feb 29;319(5867):1244-7</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">18309083</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Submicrosc Cytol Pathol. 1998 Jan;30(1):45-53</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">9530851</ArticleId></ArticleIdList></Reference><Reference><Citation>Semin ' \
                        b'Immunopathol. 2011 Sep;33(5):419-40</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">21174094</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Proteomics. 2009 Nov;9(21):4997-5000</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19810033</ArticleId></ArticleIdList></Reference><Reference><Citation>Eur J ' \
                        b'Immunol. 2001 Oct;31(10):2892-900</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">11592064</ArticleId></ArticleIdList></Reference><Reference><Citation>J Exp ' \
                        b'Med. 2002 May 20;195(10):1303-16</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">12021310</ArticleId></ArticleIdList></Reference><Reference><Citation>J Biol ' \
                        b'Chem. 1987 Jul 5;262(19):9412-20</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">3597417</ArticleId></ArticleIdList></Reference><Reference><Citation>J Biol ' \
                        b'Chem. 2009 Dec 4;284(49):34211-22</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19801663</ArticleId></ArticleIdList></Reference><Reference><Citation>Cancer ' \
                        b'Immunol Immunother. 2006 Jul;55(7):808-18</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16283305</ArticleId></ArticleIdList></Reference><Reference><Citation>J Cell ' \
                        b'Biol. 1985 Sep;101(3):942-8</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">2993317</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Proteomics Clin Appl. 2007 Nov;1(11):1446-61</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">21136642</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Transl Med. 2008 Jul 22;6:37</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">18644158</ArticleId></ArticleIdList></Reference><Reference><Citation>Nat ' \
                        b'Cell Biol. 2007 Jun;9(6):654-9</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">17486113</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Biochim Biophys Acta. 1985 Sep 9;822(2):203-18</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">2992593</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Leukemia. 2006 May;20(5):847-56</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16453000</ArticleId></ArticleIdList></Reference><Reference><Citation>Nat ' \
                        b'Rev Immunol. 2009 Aug;9(8):581-93</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19498381</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Science. 2005 Jun 24;308(5730):1862-3</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">15976285</ArticleId></ArticleIdList></Reference><Reference><Citation>Mol ' \
                        b'Cell Neurosci. 2011 Feb;46(2):409-18</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">21111824</ArticleId></ArticleIdList></Reference><Reference><Citation>J Exp ' \
                        b'Med. 1996 Mar 1;183(3):1161-72</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">8642258</ArticleId></ArticleIdList></Reference><Reference><Citation>J Cell ' \
                        b'Biol. 1999 Nov 1;147(3):599-610</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">10545503</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Biochem J. 1991 Mar 1;274 ( Pt 2):381-6</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">1848755</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Immunol. 2001 Jun 15;166(12):7309-18</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">11390481</ArticleId></ArticleIdList></Reference></ReferenceList' \
                        b'></PubmedData></PubmedArticle><PubmedArticle><MedlineCitation Status="MEDLINE" ' \
                        b'Owner="NLM"><PMID Version="1">21876761</PMID><DateCompleted><Year>2011</Year><Month>12' \
                        b'</Month><Day>30</Day></DateCompleted><DateRevised><Year>2021</Year><Month>10</Month><Day>20' \
                        b'</Day></DateRevised><Article PubModel="Print-Electronic"><Journal><ISSN ' \
                        b'IssnType="Electronic">1932-6203</ISSN><JournalIssue ' \
                        b'CitedMedium="Internet"><Volume>6</Volume><Issue>8</Issue><PubDate><Year>2011</Year' \
                        b'></PubDate></JournalIssue><Title>PloS one</Title><ISOAbbreviation>PLoS ' \
                        b'One</ISOAbbreviation></Journal><ArticleTitle>A signature in HIV-1 envelope leader peptide ' \
                        b'associated with transition from acute to chronic infection impacts envelope processing and ' \
                        b'infectivity.</ArticleTitle><Pagination><MedlinePgn>e23673</MedlinePgn></Pagination' \
                        b'><ELocationID EIdType="doi" ' \
                        b'ValidYN="Y">10.1371/journal.pone.0023673</ELocationID><Abstract><AbstractText>Mucosal ' \
                        b'transmission of the human immunodeficiency virus (HIV) results in a bottleneck in viral ' \
                        b'genetic diversity. Gnanakaran and colleagues used a computational strategy to identify ' \
                        b'signature amino acids at particular positions in Envelope that were associated either with ' \
                        b'transmitted sequences sampled very early in infection, or sequences sampled during chronic ' \
                        b'infection. Among the strongest signatures observed was an enrichment for the stable ' \
                        b'presence of histidine at position 12 at transmission and in early infection, ' \
                        b'and a recurrent loss of histidine at position 12 in chronic infection. This amino acid lies ' \
                        b'within the leader peptide of Envelope, a region of the protein that has been shown to ' \
                        b'influence envelope glycoprotein expression and virion infectivity. We show a strong ' \
                        b'association between a positively charged amino acid like histidine at position 12 in ' \
                        b'transmitted/founder viruses with more efficient trafficking of the nascent envelope ' \
                        b'polypeptide to the endoplasmic reticulum and higher steady-state glycoprotein expression ' \
                        b'compared to viruses that have a non-basic position 12 residue, a substitution that was ' \
                        b'enriched among viruses sampled from chronically infected individuals. When expressed in the ' \
                        b'context of other viral proteins, transmitted envelopes with a basic amino acid position 12 ' \
                        b'were incorporated at higher density into the virus and exhibited higher infectious titers ' \
                        b'than did non-signature envelopes. These results support the potential utility of using a ' \
                        b'computational approach to examine large viral sequence data sets for functional signatures ' \
                        b'and indicate the importance of Envelope expression levels for efficient HIV ' \
                        b'transmission.</AbstractText></Abstract><AuthorList CompleteYN="Y"><Author ' \
                        b'ValidYN="Y"><LastName>Asmal</LastName><ForeName>Mohammed</ForeName><Initials>M</Initials' \
                        b'><AffiliationInfo><Affiliation>Division of Viral Pathogenesis, Beth Israel Deaconess ' \
                        b'Medical Center, Harvard Medical School, Boston, Massachusetts, United States of America. ' \
                        b'masmal@bidmc.harvard.edu</Affiliation></AffiliationInfo></Author><Author ' \
                        b'ValidYN="Y"><LastName>Hellmann</LastName><ForeName>Ina</ForeName><Initials>I</Initials' \
                        b'></Author><Author ValidYN="Y"><LastName>Liu</LastName><ForeName>Weimin</ForeName><Initials' \
                        b'>W</Initials></Author><Author ValidYN="Y"><LastName>Keele</LastName><ForeName>Brandon ' \
                        b'F</ForeName><Initials>BF</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Perelson</LastName><ForeName>Alan ' \
                        b'S</ForeName><Initials>AS</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Bhattacharya</LastName><ForeName>Tanmoy</ForeName><Initials>T' \
                        b'</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Gnanakaran</LastName><ForeName>S</ForeName><Initials>S</Initials' \
                        b'></Author><Author ValidYN="Y"><LastName>Daniels</LastName><ForeName>Marcus</ForeName' \
                        b'><Initials>M</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Haynes</LastName><ForeName>Barton ' \
                        b'F</ForeName><Initials>BF</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Korber</LastName><ForeName>Bette ' \
                        b'T</ForeName><Initials>BT</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Hahn</LastName><ForeName>Beatrice ' \
                        b'H</ForeName><Initials>BH</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Shaw</LastName><ForeName>George ' \
                        b'M</ForeName><Initials>GM</Initials></Author><Author ' \
                        b'ValidYN="Y"><LastName>Letvin</LastName><ForeName>Norman ' \
                        b'L</ForeName><Initials>NL</Initials></Author></AuthorList><Language>eng</Language><GrantList ' \
                        b'CompleteYN="Y"><Grant><GrantID>R01 AI028433</GrantID><Acronym>AI</Acronym><Agency>NIAID NIH ' \
                        b'HHS</Agency><Country>United States</Country></Grant><Grant><GrantID>R01 ' \
                        b'RR006555</GrantID><Acronym>RR</Acronym><Agency>NCRR NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>AI28433-19</GrantID><Acronym>AI</Acronym><Agency' \
                        b'>NIAID NIH HHS</Agency><Country>United States</Country></Grant><Grant><GrantID>U19 ' \
                        b'AI067854</GrantID><Acronym>AI</Acronym><Agency>NIAID NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>P30 ' \
                        b'AI060354</GrantID><Acronym>AI</Acronym><Agency>NIAID NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>T32 ' \
                        b'AI007387</GrantID><Acronym>AI</Acronym><Agency>NIAID NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>R01 ' \
                        b'OD011095</GrantID><Acronym>OD</Acronym><Agency>NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>R01 ' \
                        b'AI094604</GrantID><Acronym>AI</Acronym><Agency>NIAID NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>RR06555-18</GrantID><Acronym>RR</Acronym><Agency' \
                        b'>NCRR NIH HHS</Agency><Country>United States</Country></Grant><Grant><GrantID>U01 ' \
                        b'AI067854</GrantID><Acronym>AI</Acronym><Agency>NIAID NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>R37 ' \
                        b'AI028433</GrantID><Acronym>AI</Acronym><Agency>NIAID NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant><Grant><GrantID>AI-067854</GrantID><Acronym>AI</Acronym><Agency' \
                        b'>NIAID NIH HHS</Agency><Country>United ' \
                        b'States</Country></Grant></GrantList><PublicationTypeList><PublicationType ' \
                        b'UI="D016428">Journal Article</PublicationType><PublicationType UI="D052061">Research ' \
                        b'Support, N.I.H., Extramural</PublicationType><PublicationType UI="D013485">Research ' \
                        b'Support, Non-U.S. Gov\'t</PublicationType></PublicationTypeList><ArticleDate ' \
                        b'DateType="Electronic"><Year>2011</Year><Month>08</Month><Day>18</Day></ArticleDate' \
                        b'></Article><MedlineJournalInfo><Country>United States</Country><MedlineTA>PLoS ' \
                        b'One</MedlineTA><NlmUniqueID>101285081</NlmUniqueID><ISSNLinking>1932-6203</ISSNLinking' \
                        b'></MedlineJournalInfo><ChemicalList><Chemical><RegistryNumber>0</RegistryNumber' \
                        b'><NameOfSubstance UI="D015686">Gene Products, ' \
                        b'env</NameOfSubstance></Chemical><Chemical><RegistryNumber>0</RegistryNumber' \
                        b'><NameOfSubstance UI="D016655">HIV Core Protein ' \
                        b'p24</NameOfSubstance></Chemical><Chemical><RegistryNumber>0</RegistryNumber' \
                        b'><NameOfSubstance UI="D021382">Protein Sorting ' \
                        b'Signals</NameOfSubstance></Chemical></ChemicalList><CitationSubset>IM</CitationSubset' \
                        b'><MeshHeadingList><MeshHeading><DescriptorName UI="D000208" MajorTopicYN="N">Acute ' \
                        b'Disease</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D000595" ' \
                        b'MajorTopicYN="N">Amino Acid ' \
                        b'Sequence</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D002908" ' \
                        b'MajorTopicYN="N">Chronic Disease</DescriptorName></MeshHeading><MeshHeading><DescriptorName ' \
                        b'UI="D018450" MajorTopicYN="Y">Disease ' \
                        b'Progression</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D004721" ' \
                        b'MajorTopicYN="N">Endoplasmic Reticulum</DescriptorName><QualifierName UI="Q000378" ' \
                        b'MajorTopicYN="N">metabolism</QualifierName></MeshHeading><MeshHeading><DescriptorName ' \
                        b'UI="D015686" MajorTopicYN="N">Gene Products, env</DescriptorName><QualifierName ' \
                        b'UI="Q000737" MajorTopicYN="N">chemistry</QualifierName><QualifierName UI="Q000378" ' \
                        b'MajorTopicYN="Y">metabolism</QualifierName></MeshHeading><MeshHeading><DescriptorName ' \
                        b'UI="D016655" MajorTopicYN="N">HIV Core Protein p24</DescriptorName><QualifierName ' \
                        b'UI="Q000378" MajorTopicYN="N">metabolism</QualifierName></MeshHeading><MeshHeading' \
                        b'><DescriptorName UI="D015658" MajorTopicYN="N">HIV ' \
                        b'Infections</DescriptorName><QualifierName UI="Q000473" ' \
                        b'MajorTopicYN="Y">pathology</QualifierName><QualifierName UI="Q000821" ' \
                        b'MajorTopicYN="Y">virology</QualifierName></MeshHeading><MeshHeading><DescriptorName ' \
                        b'UI="D015497" MajorTopicYN="N">HIV-1</DescriptorName><QualifierName UI="Q000378" ' \
                        b'MajorTopicYN="Y">metabolism</QualifierName><QualifierName UI="Q000472" ' \
                        b'MajorTopicYN="Y">pathogenicity</QualifierName></MeshHeading><MeshHeading><DescriptorName ' \
                        b'UI="D006801" MajorTopicYN="N">Humans</DescriptorName></MeshHeading><MeshHeading' \
                        b'><DescriptorName UI="D019169" MajorTopicYN="N">Jurkat ' \
                        b'Cells</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D016013" ' \
                        b'MajorTopicYN="N">Likelihood ' \
                        b'Functions</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D008954" ' \
                        b'MajorTopicYN="N">Models, ' \
                        b'Biological</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D008969" ' \
                        b'MajorTopicYN="N">Molecular Sequence ' \
                        b'Data</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D010802" ' \
                        b'MajorTopicYN="N">Phylogeny</DescriptorName></MeshHeading><MeshHeading><DescriptorName ' \
                        b'UI="D011110" MajorTopicYN="N">Polymorphism, ' \
                        b'Genetic</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D021382" ' \
                        b'MajorTopicYN="Y">Protein Sorting ' \
                        b'Signals</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D017434" ' \
                        b'MajorTopicYN="N">Protein Structure, ' \
                        b'Tertiary</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D021381" ' \
                        b'MajorTopicYN="N">Protein ' \
                        b'Transport</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D014771" ' \
                        b'MajorTopicYN="N">Virion</DescriptorName><QualifierName UI="Q000378" ' \
                        b'MajorTopicYN="N">metabolism</QualifierName></MeshHeading></MeshHeadingList' \
                        b'></MedlineCitation><PubmedData><History><PubMedPubDate ' \
                        b'PubStatus="received"><Year>2011</Year><Month>05</Month><Day>24</Day></PubMedPubDate' \
                        b'><PubMedPubDate PubStatus="accepted"><Year>2011</Year><Month>07</Month><Day>22</Day' \
                        b'></PubMedPubDate><PubMedPubDate ' \
                        b'PubStatus="entrez"><Year>2011</Year><Month>8</Month><Day>31</Day><Hour>6</Hour><Minute>0' \
                        b'</Minute></PubMedPubDate><PubMedPubDate ' \
                        b'PubStatus="pubmed"><Year>2011</Year><Month>8</Month><Day>31</Day><Hour>6</Hour><Minute>0' \
                        b'</Minute></PubMedPubDate><PubMedPubDate ' \
                        b'PubStatus="medline"><Year>2011</Year><Month>12</Month><Day>31</Day><Hour>6</Hour><Minute>0' \
                        b'</Minute></PubMedPubDate></History><PublicationStatus>ppublish</PublicationStatus' \
                        b'><ArticleIdList><ArticleId IdType="pubmed">21876761</ArticleId><ArticleId ' \
                        b'IdType="doi">10.1371/journal.pone.0023673</ArticleId><ArticleId ' \
                        b'IdType="pii">PONE-D-11-09221</ArticleId><ArticleId ' \
                        b'IdType="pmc">PMC3158090</ArticleId></ArticleIdList><ReferenceList><Reference><Citation>J ' \
                        b'Virol. 2002 Jun;76(11):5315-25</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">11991960</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Virology. 1992 Jul;189(1):103-10</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">1376536</ArticleId></ArticleIdList></Reference><Reference><Citation>AIDS. ' \
                        b'2003 Sep 5;17(13):1871-9</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">12960819</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 1998 Apr;72(4):2855-64</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">9525605</ArticleId></ArticleIdList></Reference><Reference><Citation>Mol ' \
                        b'Med. 2006 Jul-Aug;12(7-8):137-42</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">17088945</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Science. 1996 Mar 15;271(5255):1582-6</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">8599114</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 2007 Oct;81(20):10869-78</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">17670815</ArticleId></ArticleIdList></Reference><Reference><Citation>Mol ' \
                        b'Biol Cell. 2005 Jan;16(1):279-91</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">15496459</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 2009 Apr;83(8):3556-67</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19193811</ArticleId></ArticleIdList></Reference><Reference><Citation>Mol ' \
                        b'Biol Evol. 2010 Feb;27(2):221-4</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19854763</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 2011 Jul 12;108(28):11530-5</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">21690407</ArticleId></ArticleIdList></Reference><Reference><Citation>Int ' \
                        b'Rev Cytol. 2005;245:91-121</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16125546</ArticleId></ArticleIdList></Reference><Reference><Citation>FEBS ' \
                        b'Lett. 2006 Jun 26;580(15):3775-8</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16777098</ArticleId></ArticleIdList></Reference><Reference><Citation>Nature' \
                        b'. 1995 Jan 12;373(6510):123-6</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">7816094</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 1993 Jun;67(6):3345-56</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">8497055</ArticleId></ArticleIdList></Reference><Reference><Citation>Trends ' \
                        b'Biochem Sci. 2006 Oct;31(10):563-71</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16919958</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 1992 Nov 1;89(21):10247-51</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">1438212</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Virology. 1994 Oct;204(1):266-78</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">8091657</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 1996 Sep 3;93(18):9606-11</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">8790377</ArticleId></ArticleIdList></Reference><Reference><Citation>Science' \
                        b'. 2004 Mar 26;303(5666):2019-22</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">15044802</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 2002 Dec;76(23):11953-9</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">12414937</ArticleId></ArticleIdList></Reference><Reference><Citation>J Cell ' \
                        b'Sci. 2006 Nov 1;119(Pt 21):4373-80</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">17074831</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Virology. 2000 Jul 5;272(2):417-28</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">10873786</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Theor Biol. 2000 Apr 7;203(3):285-301</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">10716909</ArticleId></ArticleIdList></Reference><Reference><Citation>Cell. ' \
                        b'2006 Dec 1;127(5):999-1013</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">17129784</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Science. 2007 Mar 16;315(5818):1583-6</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">17363674</ArticleId></ArticleIdList></Reference><Reference><Citation>FASEB ' \
                        b'J. 2003 Jun;17(9):1058-67</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">12773488</ArticleId></ArticleIdList></Reference><Reference><Citation>BMC ' \
                        b'Evol Biol. 2006 Mar 23;6:28</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16556318</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 2003 Dec 23;100(26):15812-7</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">14668432</ArticleId></ArticleIdList></Reference><Reference><Citation>J Exp ' \
                        b'Med. 2009 Jun 8;206(6):1253-72</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">19487423</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 1997 Oct;71(10):7518-25</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">9311831</ArticleId></ArticleIdList></Reference><Reference><Citation>Nat Rev ' \
                        b'Microbiol. 2006 Apr;4(4):312-7</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16541138</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 2006 Sep 19;103(38):13950-5</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16966601</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 1990 Dec;64(12):6297-304</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">2243395</ArticleId></ArticleIdList></Reference><Reference><Citation>Nature. ' \
                        b'2003 Mar 20;422(6929):307-12</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">12646921</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Virology. 2005 Feb 5;332(1):418-29</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">15661172</ArticleId></ArticleIdList></Reference><Reference><Citation>Curr ' \
                        b'Top Microbiol Immunol. 2005;285:175-98</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">15609504</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 1992 Apr;66(4):2296-301</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">1548763</ArticleId></ArticleIdList></Reference><Reference><Citation>Arch ' \
                        b'Virol. 1997;142(1):37-51</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">9155871</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 2006 Jul;80(14):7226-34</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">16809328</ArticleId></ArticleIdList></Reference><Reference><Citation>Proc ' \
                        b'Natl Acad Sci U S A. 2008 May 27;105(21):7552-7</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">18490657</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 1999 Dec;73(12):10489-502</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">10559367</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Science. 1992 Feb 28;255(5048):1134-7</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">1546316</ArticleId></ArticleIdList></Reference><Reference><Citation>Protein ' \
                        b'Expr Purif. 2008 Oct;61(2):142-8</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">18595733</ArticleId></ArticleIdList></Reference><Reference><Citation>Endocr ' \
                        b'Rev. 2008 May;29(3):317-33</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">18436705</ArticleId></ArticleIdList></Reference><Reference><Citation' \
                        b'>Science. 1993 Aug 27;261(5125):1179-81</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">8356453</ArticleId></ArticleIdList></Reference><Reference><Citation>J ' \
                        b'Virol. 2000 Sep;74(18):8358-67</Citation><ArticleIdList><ArticleId ' \
                        b'IdType="pubmed">10954535</ArticleId></ArticleIdList></Reference></ReferenceList' \
                        b'></PubmedData></PubmedArticle></PubmedArticleSet> '


def test_get_papers_from_NCBI():
    ncbi_m.get_papers_from_NCBI("mirna")


def test_efetch_manula():
    """
    This test is to check the manual made efetch that is a copy with a sligthly modification to the
    one in eutilz.
    :return:
    """
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
    eutilcito.get_ids_information(db_id=[35521437, 35524416, 21876726], db='pubmed')


def test_parse_pubmed_data():
    eutilcito.get_ids_information(db_id=[35521437, 35524416, 21876726], db='pubmed')


def test_parse_pubmed_no_year():
    eutilcito.get_ids_information(db_id=[35528540, 34662018], db='pubmed')


def test_parse_pubmed_data_book_tons_cites():
    egs = book_info_xml.decode('latin')
    egs = ET.fromstring(egs)
    cites = tons_cites_xml.decode('latin')
    cites = ET.fromstring(cites)
    r = eutilcito.parse_pubmed_data(egs=egs, cites=cites, pid='123')
    article_data = {'Pubmed_id': '123',
                    'Year': '2022',
                    'Title': 'StatPearls',
                    'Key_words': [],
                    'Cite_by': 373}
    assert r == article_data, f'The book had xml with the values {article_data}, but {r} was retrived.'


def test_parse_pubmed_data_book_no_cites():
    egs = book_info_one_xml.decode('latin')
    egs = ET.fromstring(egs)
    cites = no_cites_one_xml.decode('latin')
    cites = ET.fromstring(cites)
    r = eutilcito.parse_pubmed_data(egs=egs, cites=cites, pid='123')
    article_data = {'Pubmed_id': '123',
                    'Year': '2022',
                    'Title': 'StatPearls',
                    'Key_words': [],
                    'Cite_by': 0}
    assert r == article_data, f'The book had xml with the values {article_data}, but {r} was retrived.'


def test_multiple_cites():
    client = Client(api_key='8fa634394aad95b85063646ab64bef2c2408')
    r = eut._multiple_elink(self=client, db='pubmed', ids=[21876726, 21876761])
    print(r)


def test_multiple_fetch_info():
    client = Client(api_key='8fa634394aad95b85063646ab64bef2c2408')
    r = eut._multiple_efetch(self=client, db='pubmed', ids=[21876726, 21876761])
    cites = no_cites_one_xml.decode('latin')
    cites = ET.fromstring(cites)


def test_multiple_fetch_info():
    egs = article_info_many_xml.decode('latin')
    egs = ET.fromstring(egs)
    r = eutilcito.parse_many_pubmed_data(egs=egs)
    article_data = [{'Pubmed_id': '21876726',
                     'Year': '2011',
                     'Title': 'Exosomes: secreted vesicles and intercellular communications.',
                     'Key_words': []},
                    {'Pubmed_id': '21876761',
                     'Year': '2011',
                     'Title': 'A signature in HIV-1 envelope leader peptide associated with transition from acute to '
                              'chronic infection impacts envelope processing and infectivity.',
                     'Key_words': []}
                    ]
    assert r == article_data, f"No"


def test_multiple_cites():
    cites = ET.fromstring(cites_many_xml.decode('latin'))
    r = eutilcito.parse_many_cites(cites=cites)
    article_data = [{'Pubmed_id': '21876726',
                     'Cite_by': 374},
                    {'Pubmed_id': '21876761',
                     'Cite_by': 43}
                    ]
    assert r == article_data, f"Cite data don't match"
    print(r)


def test_join_info_w_cites():
    cites = [{'Pubmed_id': '21876726',
              'Cite_by': 374},
             {'Pubmed_id': '21876761',
              'Cite_by': 43}
             ]
    info = [{'Pubmed_id': '21876726',
             'Year': '2011',
             'Title': 'Exosomes: secreted vesicles and intercellular communications.',
             'Key_words': []},
            {'Pubmed_id': '21876761',
             'Year': '2011',
             'Title': 'A signature in HIV-1 envelope leader peptide associated with transition from acute to '
                      'chronic infection impacts envelope processing and infectivity.',
             'Key_words': []}
            ]
    x = eutilcito.join_info_w_cites(info=info, cites=cites)
    print(x)


def test_get_papers_information():
    information = eutilcito.get_papers_information([21876726, 21876761])
    guide = [{'Pubmed_id': '21876726', 'Year': '2011', 'Title': 'Exosomes: secreted vesicles and intercellular '
                                                                'communications.', 'Key_words': [], 'Cite_by': 374},
             {'Pubmed_id': '21876761', 'Year': '2011', 'Title': 'A signature in HIV-1 envelope leader peptide '
                                                                'associated with transition from acute to chronic '
                                                                'infection impacts envelope processing and '
                                                                'infectivity.', 'Key_words': [], 'Cite_by': 43}]
    assert information == guide, f"Something bad"
