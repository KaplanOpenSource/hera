# -*- coding: utf-8 -*-
"""
    Converts the bib files to hebrew.
"""
import unicodedata

import os


class bibItem:
    """
        Represents a single bib item.

        \bibitem{paper2}
            John Smith and Keren Wilberforth.
        \newblock Is $4\frac{m}{s}$ the right wind speed.
        \newblock {\em Journal of internal affairs}, October 2011.
    """

    _item = None

    def __init__(self,bibItem):
        """
            Initializes the object.
        :param bibItem: list
            A list of strings.
        """
        self._item = bibItem

    @property
    def isHebrew(self):
        """
            True if a hebrew character is detected in the _item .
        :return: bool.
        """
        return any([self.isHebrewText(line) for line in self._item])


    def isHebrewText(self,line):
        return any(['HEBREW' in self._getUnicodeName(c) for c in line])


    def _getUnicodeName(self,c):
        try:
            return unicodedata.name(c)
        except ValueError:
            return 'NONE'

    def convert(self):
        if self.isHebrew:
            newItem = ['\\sethebrew\n']
            for line in self._item:
                # Now we need to scan each word:
                #   * add \L{} around english words.
                #   * Replace AND with 'ן'.
                #   * add $..$ around numbers.

                # 1. split to words.
                sentence = []
                addVe = False
                for indx,word in enumerate(line.split(" ")):
                    if word.startswith("\\"):
                        sentence.append(word)
                        continue

                    if word.strip().lower() == "and":
                        # if there is an item before this one,
                        # and it has comma in the end, remove it.
                        # raise a flag to and the 've' to the next word.
                        if len(sentence) > indx-1:
                            if sentence[indx-1][-1] ==',':
                                # remove the comma.
                                sentence[indx - 1] = sentence[indx - 1][:-1]
                            addVe = True
                            continue

                    theword = [] if not addVe else ['ו']

                    indx = 0
                    while indx < len(word):
                        charName = self._getUnicodeName(word[indx])

                        if 'LATIN' in charName:
                            # begin an english word, now search for its end.
                            beginIndex = indx
                            theword.append('\\L{')
                            while 'LATIN' in self._getUnicodeName(word[indx]):
                                indx += 1
                                if indx == len(word):
                                    break

                            theword.append(word[beginIndex:indx])
                            theword.append('}')

                        elif  word[indx] == '\\':
                            # begin command:
                            beginIndex = indx
                            indx += 1
                            if indx == len(word):
                                theword.append("\\")
                            else:
                                while 'LATIN' in self._getUnicodeName(word[indx]):
                                    indx += 1
                                    if indx == len(word):
                                        break
                                theword.append(word[beginIndex:indx])
                        elif  word[indx] == 'DOLLAR SIGN':
                            # begin command:
                            beginIndex = indx
                            indx += 1
                            if indx == len(word):
                                theword.append("$")
                            else:
                                while 'DOLLAR SIGN' in self._getUnicodeName(word[indx]):
                                    indx += 1
                                    if indx == len(word):
                                        break
                                theword.append(word[beginIndex:indx])


                        elif 'DIGIT'  in charName:
                            # begin number
                            beginIndex = indx
                            theword.append('$')
                            charType = self._getUnicodeName(word[indx])
                            while ('DIGIT' in charType) or ('SOLIDUS' in charType):
                                indx += 1
                                if indx == len(word):
                                    break
                                charType = self._getUnicodeName(word[indx])

                            theword.append(word[beginIndex:indx])
                            theword.append('$')
                        else:
                            theword.append(word[indx])
                            indx +=1

                    addVe = False
                    sentence.append("".join(theword))
                newItem.append(" ".join(sentence))
        else:
            newItem = ['\\unsethebrew\n'] + self._item
        return "".join(newItem)


class bibtexFile:
    """
        Parses the bibItems and constructs them.
    """

    _bibItems = None
    _first_last_Lines = None

    def __init__(self,thefile):

        thefileStr = None

        if isinstance(thefile,str):
            if os.path.exists(thefile):
                with open(thefile,"r") as inputFile:
                    thefileStr = inputFile.readlines()
            else:
                raise ValueError(f"The file {thefile} does not exist")
        else:
            thefileStr = thefile.readlines()

        self._bibItems = []

        self._first_last_Lines = (thefileStr[0],thefileStr[-1])
        self._parseBibItems(thefileStr[1:-1])

    def _parseBibItems(self,strList):
        """
            Traverse over the list and populates the _bibItems list.
            An item is between bibitem.

        :param strList: list
                a list of bibtex file.
        :return:
                None
        """

        bibitemlines = []

        inItem = False

        curItem = []
        for line in strList:
            if line.startswith("\\bibitem"):
                #close current item
                self._bibItems.append(bibItem(curItem))
                curItem = [line]
            else:
                curItem.append(line)

        self._bibItems.append(bibItem(curItem))

    @property
    def items(self):
        return self._bibItems

    def convert(self):
        final = []
        for itm in bb.items:
            final.append(itm.convert())

        return "\n".join([self._first_last_Lines[0]] + final + [self._first_last_Lines[0]])


if __name__ == "__main__":
    bb = bibtexFile("/home/yehudaa/Projects/2020/hebrewCrossTex/output.bbl")
    print(bb.convert())



officcrop
