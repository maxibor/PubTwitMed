# PubTwitMed

[![Build Status](https://travis-ci.org/maxibor/PubTwitMed.svg?branch=master)](https://travis-ci.org/maxibor/PubTwitMed)

#### A twitter bot publishing regularly all recent scientific articles from PubMed related to a given subject.


### Requirements
- Python 3.x
- [BioPython](https://biopython.org/wiki/Biopython)
- [Tweepy](http://www.tweepy.org/)

--

### Manual

To get twitter API credentials, please have a look at the documentation of [Tweepy](http://docs.tweepy.org/en/latest/getting_started.html).

#### PubTwitMed Usage :
```
$ python pubmed_twitter_bot.py -h
usage: pubtwitmed [-h] [-doi DOI_DB] [-artmax ARTMAX] [-topic TOPIC]
                  [-email ENTREZ_EMAIL] [-ck CONSUMER_KEY]
                  [-cs CONSUMER_SECRET] [-at ACCESS_TOKEN]
                  [-ats ACCESS_TOKEN_SECRET] [-nak NCBI_API_KEY]

==========================================================
PubTwitMed
A twitter bot publishing regularly all recent scientific
articles from PubMed related to a given subject.
Author: Maxime Borry
Contact: <maxime.borry[at]gmail.com>
Homepage & Documentation: github.com/maxibor/PubTwitMed
==========================================================


optional arguments:
  -h, --help            show this help message and exit
  -doi DOI_DB           Path to DOI database file /path/to/doi_db.txt
  -artmax ARTMAX        Max articles to retrieve at each query. Default = 15
  -topic TOPIC          Topic to tweet about
  -email ENTREZ_EMAIL   Email for NCBI Entrez.
  -ck CONSUMER_KEY      Twitter consumer key
  -cs CONSUMER_SECRET   Twitter consumer secret
  -at ACCESS_TOKEN      Twitter access token
  -ats ACCESS_TOKEN_SECRET
                        Twitter access token secret
  -nak NCBI_API_KEY     NCBI api key
```


- *To execute it regularly, please have a look at [Crontab](https://en.wikipedia.org/wiki/Cron)*
