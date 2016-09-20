#!/usr/bin/env python
import time

topic_to_tweet_about = "microbiome"
max_number_of_articles_to_retrieve = 15

def twitterbot(string_to_tweet) :
    import tweepy
    login = tweepy.OAuthHandler("consumer_key", "consumer_secret")
    login.set_access_token("access_token", "access_token_secret")
    this_api = tweepy.API(login)
    this_api.update_status(status=string_to_tweet)

def pubmed_search(search_term,nb_max_articles) :
    '''
    Search Pubmed for the nb_max_articles most recent articles on the
    search_term subject.

    INPUT : Search Term(str) and nb_max_articles(int)
    OUPUT : Dictionnary of Lists ['DOI':['Title','First Author','PubDate']]
    EXAMPLE : pubmed_search("microbiome",10)
    '''
    from Bio import Entrez

    article_dictionary = {}
    Entrez.email="maxime.borry@gmail.com"
    max_number_of_articles = nb_max_articles

    myhandle = Entrez.esearch(db="pubmed",term=search_term,
    retmax = max_number_of_articles)
    my_record = Entrez.read(myhandle)
    my_list = my_record["IdList"]
    for article in range(0, len(my_list)):
        listId = my_list[article]
        my_secondary_handle = Entrez.esummary(db="pubmed", id=listId)
        my_record = Entrez.read(my_secondary_handle)
        # for akey in my_record[0].keys() : #list of keys
        #     print akey, ":", my_record[0][akey]
        # print my_record[0].keys()
        #need "DOI"
        #need "Title"
        #need "AuthorList[0]"
        #need "PubDate"
        one_article = my_record[0]
        try :
            article_dictionary[one_article["DOI"]] = [one_article["Title"],
            one_article["AuthorList"][0],one_article["PubDate"]]
        except :
            continue
        time.sleep(0.5)
        # break
    return(article_dictionary)

def doi_tool (adoi) :
    '''
    Gets a DOI in input, and check if it already in doi_db.txt file.
    If yes returns "already", if not appends it to doi_db.txt and returns the
    direct link to the publication

    INPUT : DOI_identifier (string)
    OUTPUT : DOI url
    EXAMPLE : doi_tool ("10.1371/journal.pone.0161211")
    '''

    def doi_url_resolver (adoi) :
        '''
        Gets a DOI in input and uses the dx.doi.org service to
        return article url
        '''
        return("http://dx.doi.org/"+str(adoi))

    def doi_checker(doi_string) :
        '''
        Gets a DOI in input, and check if it already in doi_db.txt file.
        If yes, returns TRUE, if not, returns FALSE and DOI
        '''
        dois_list=[]
        with open("doi_db.txt","r") as doi_db :
            for line in doi_db :
                line=line.rstrip()
                dois_list.append(line)
        if doi_string in dois_list :
            return(True,"NA")
        elif doi_string not in dois_list :
            dois_list.append(doi_string)
            return(False,doi_string)

    def doi_appender(doi_string) :
        '''
        Appends DOI in doi_db.txt
        '''
        with open("doi_db.txt","a") as doi_db :
                doi_db.write(doi_string+"\n")



    if doi_checker(adoi)[0] == False :
        doi_appender(doi_checker(adoi)[1])
        return(doi_url_resolver(adoi))
    else :
        return("already")

def url_shortener (full_url) :
    from pyshorteners import Shortener

    url = full_url
    shortener = Shortener('Isgd')
    return(shortener.short(url))

def string_shortener (string_to_shorten,max_size) :
    '''
    Shortens titles strings that are more than max_size
    Returns shortened titled strings

    INPUT : title_string
    OUPUT : str(title_string)+"..."
    EXAMPLE : title_shortener(title)
    '''
    return string_to_shorten[0:max_size]+"..."

myquery =  pubmed_search(topic_to_tweet_about,max_number_of_articles_to_retrieve)
for article in myquery :
    mystatus = doi_tool(article)
    if mystatus != "already" :
        print "DOI : ", article
        print "URL : ", mystatus
        print "Title : ", myquery[article][0]
        # final_title = string_shortener(myquery[article][0],60)
        print "First Author : ", myquery[article][1]
        final_author = myquery[article][1]+" et al."
        print "PubDate : ", myquery[article][2]
        final_date = myquery[article][2]
        hashtag = "#microbiome"
        try :
            print "Short-URL : ", url_shortener(mystatus)
            final_url = url_shortener(mystatus)
        except :
            continue
        try :
            almost_to_tweet =" - "+final_author+" - "+final_date+" - "+final_url+" "+hashtag
            max_title_len = 130 - len(almost_to_tweet) #137 because of "..." after title
            final_title = string_shortener(myquery[article][0],max_title_len)
            text_to_tweet = final_title+" - "+final_author+" - "+final_date+" - "+final_url+" "+hashtag
            print text_to_tweet
            print len(text_to_tweet)
            twitterbot(text_to_tweet)
            sleep(10)

        except :
            continue
    break
