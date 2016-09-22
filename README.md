#PubTwitMed
####A twitter bot publishing regularly all recent scientific articles from PubMed related to a given subject.


###Requirements
- Python 2.7.x
- [BioPython](https://biopython.org/wiki/Biopython)
- [Tweepy](http://www.tweepy.org/)
- [Pyshorteners](https://github.com/ellisonleao/pyshorteners)

--

###Manual
- The article subject to tweet about is to be edited in the variable **topic_to_tweet_about**

- To add you twitter API credentials, please have a look at the documentation of [Tweepy](http://tweepy.readthedocs.io/en/v3.5.0/getting_started.html#hello-tweepy)

- #####PubTwitMed Usage :
`./pubmed_twitter_bot.py`

- *To execute it regularly, please have a look at [Crontab](https://en.wikipedia.org/wiki/Cron)*

	Contrab file example :
	`@daily /path/to/directory/pubmed_twitter_bot.py >> path/to/logfile.txt 2>&1`
