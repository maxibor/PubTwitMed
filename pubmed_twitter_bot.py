#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pubmed_twitter_bot.py : a twitter bot searching for the latest articles
on a given subjet on PubMed and publishing on twitter."""

__author__ = "Maxime Borry"
__license__ = "BEERWARE"

import time
import datetime
import argparse


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='pubtwitmed',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f'''
==========================================================
PubTwitMed
A twitter bot publishing regularly all recent scientific
articles from PubMed related to a given subject.
Author: Maxime Borry
Contact: <maxime.borry[at]gmail.com>
Homepage & Documentation: github.com/maxibor/PubTwitMed
==========================================================
        ''')
    parser.add_argument(
        '-doi',
        dest="doi_db",
        default="./doi_db.txt",
        help="Path to DOI database file /path/to/doi_db.txt")
    parser.add_argument(
        '-artmax',
        dest="artmax",
        default=15,
        help="Max articles to retrieve at each query. Default = 15")
    parser.add_argument(
        '-topic',
        default=None,
        help="Topic to tweet about")
    parser.add_argument(
        '-email',
        dest="entrez_email",
        default=None,
        help="Email for NCBI Entrez.")
    parser.add_argument(
        '-ck',
        dest="consumer_key",
        default=None,
        help="Twitter consumer key"
    )
    parser.add_argument(
        '-cs',
        dest="consumer_secret",
        default=None,
        help="Twitter consumer secret"
    )
    parser.add_argument(
        '-at',
        dest="access_token",
        default=None,
        help="Twitter access token"
    )
    parser.add_argument(
        '-ats',
        dest="access_token_secret",
        default=None,
        help="Twitter access token secret"
    )

    args = parser.parse_args()

    doi_db = args.doi_db
    art_max = int(args.artmax)
    topic = args.topic
    entrez_email = args.entrez_email
    cons_key = args.consumer_key
    cons_secret = args.consumer_secret
    acc_tok = args.access_token
    acc_tok_sec = args.access_token_secret

    return(doi_db, art_max, topic, entrez_email, cons_key, cons_secret, acc_tok, acc_tok_sec)


TOPIC_TO_SEARCH_AND_TWEET = "a topic to search about on pubmed"
MAX_NB_ART_TO_GET = 15
PATH_TO_DOI_DB = "/path/to/doi_db.txt"


def twitterbot(string_to_tweet, ck, cs, at, ats):
    '''
    Publish on twitter the string_to_tweet

    INPUT :
        string_to_tweet(str): text to tweet
        ck(str): Twitter consumer key
        cs(str): Twitter consumer secret
        at(str): Twitter access token
        ats(str): Twitter access token secret

    OUTPUT : None
    '''
    import tweepy
    login = tweepy.OAuthHandler(ck, cs)
    login.set_access_token(at, ats)
    this_api = tweepy.API(login)
    this_api.update_status(status=string_to_tweet)


def pubmed_search(search_term, nb_max_articles, entrez_email):
    '''
    Search Pubmed for the nb_max_articles most recent articles on the
    search_term subject.

    INPUT : Search Term(str) and nb_max_articles(int)
    OUPUT : Dictionnary of Lists ['DOI':['Title','First Author','PubDate']]
    '''
    from Bio import Entrez

    article_dictionary = {}
    Entrez.email = entrez_email
    max_number_of_articles = nb_max_articles

    myhandle = Entrez.esearch(db="pubmed", term=search_term,
                              retmax=max_number_of_articles)
    my_record = Entrez.read(myhandle)
    my_list = my_record["IdList"]
    for article in range(0, len(my_list)):
        listId = my_list[article]
        my_secondary_handle = Entrez.esummary(db="pubmed", id=listId)
        my_record = Entrez.read(my_secondary_handle)
        one_article = my_record[0]
        if len(one_article["AuthorList"]) > 1:
            authorlist = one_article["AuthorList"][0].encode(
                'utf-8').decode("utf-8")+" et al."
        else:
            authorlist = one_article["AuthorList"][0].encode(
                'utf-8').decode("utf-8")
        try:
            article_dictionary[one_article["DOI"]] = [one_article["Title"],
                                                      authorlist, one_article["PubDate"]]
        except:
            continue
        time.sleep(0.5)
        # break
    return(article_dictionary)


def string_shortener(string_to_shorten, max_size):
    '''
    Shortens titles strings that are more than max_size
    Returns shortened titled strings

    INPUT : title_string,max_size(str,int)
    OUPUT : shortened_title_string+"..."(str)
    EXAMPLE : title_shortener(title,40)
    '''
    if len(string_to_shorten) > max_size:
        return((string_to_shorten[0:max_size]+"...").capitalize())

    return(string_to_shorten.capitalize())


def doi_tool(adoi, doi_db):
    '''
    Gets a DOI in input, and check if it already in doi_db.txt file.
    If yes returns "already", if not appends it to doi_db.txt and returns the
    direct link to the publication

    INPUT :
    DOI_identifier (string)
    doi_db(str): Path to doi.txt file
    OUTPUT : DOI url
    EXAMPLE : doi_tool ("10.1371/journal.pone.0161211")
    '''

    def doi_url_resolver(adoi):
        '''
        Gets a DOI in input and uses the dx.doi.org service to
        return article url

        INPUT : DOI(str)
        OUTPUT : DOI URL(str)
        EXAMPLE : adoi("10.1371/example.doi.42")
        '''
        return("http://dx.doi.org/"+str(adoi))

    def doi_checker(doi_string, doi_db):
        '''
        Gets a DOI in input, and check if it already in doi_db.txt file.
        If yes, returns TRUE, if not, returns FALSE and DOI

        INPUT :
            DOI(str)
            doi_db(str): Path to doi.txt file
        OUTPUT : True/False(Bool)
        EXAMPLE : doi_checker("10.1371/example.doi.42")
        '''
        dois_list = []
        with open(doi_db, "r") as doi_db:
            for line in doi_db:
                line = line.rstrip()
                dois_list.append(line)
        if doi_string in dois_list:
            return(True, "NA")
        elif doi_string not in dois_list:
            dois_list.append(doi_string)
            return(False, doi_string)

    def doi_appender(doi_string, doi_db):
        '''
        Appends DOI in doi_db.txt

        INPUT :
            DOI(str)
            doi_db(str): Path to doi.txt file
        OUTPUT : None
        EXAMPLE : doi_appender("10.1371/example.doi.42")
        '''
        with open(doi_db, "a") as doi_db:
            doi_db.write(doi_string+"\n")

    if doi_checker(adoi, doi_db)[0] == False:
        doi_appender(doi_checker(adoi, doi_db)[1], doi_db)
        return(doi_url_resolver(adoi))
    else:
        return("already")


if __name__ == '__main__':
    DOI_DB, ART_MAX, TOPIC, ENTREZ_EMAIL, CONS_KEY, CONS_SECRET, ACC_TOK, ACC_TOK_SEC = _get_args()
    print("> > > > "+str(datetime.datetime.now()))
    myquery = pubmed_search(TOPIC, ART_MAX, ENTREZ_EMAIL)

    for article in myquery:
        mystatus = doi_tool(article, DOI_DB)
        if mystatus != "already":
            print("DOI : ", article)
            print("URL : ", mystatus)
            print("Title : ", myquery[article]
                  [0].encode('utf-8').decode("utf-8"))
            # final_title = string_shortener(myquery[article][0],60)
            print("First Author : ",
                  myquery[article][1])
            final_author = myquery[article][1]
            print("PubDate : ", myquery[article][2])
            final_date = myquery[article][2]
            final_url = mystatus
            hashtag = f"#{TOPIC}"
            try:
                almost_to_tweet = " - "+final_author+" - "+final_url+" "+hashtag
                max_title_len = 200 - len(almost_to_tweet)
                final_title = string_shortener(
                    myquery[article][0].encode('utf-8').decode("utf-8"), max_title_len)
                text_to_tweet = final_title+" - "+final_author+" - "+final_url+" "+hashtag
                print(text_to_tweet)
                print("tweet length :", len(text_to_tweet))
                print("= = = = = = =")
                if CONS_KEY and CONS_SECRET and ACC_TOK and ACC_TOK_SEC:
                    twitterbot(text_to_tweet, CONS_KEY,
                               CONS_SECRET, ACC_TOK, ACC_TOK_SEC)
                time.sleep(10)
            except Exception as e:
                print(e)
                continue
