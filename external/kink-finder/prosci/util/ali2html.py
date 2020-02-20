#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Sebastian Kelm
# Moved from ali.py on 07/02/2010
#

from array import array


def makeHTMLrow(title, desc, htmlseq):
    if len(title) and title[0]=='>':
      title=title[1:]
    title=title.replace('>', '&gt;').replace('<', '&lt;')
    desc=desc.replace('>', '&gt;').replace('<', '&lt;')
    return '<tr class="seqalirow"><td class="seqtitle">%s</td><td class="seqdesc">%s</td><td class="seq">%s</td></tr>' % (title, desc, htmlseq)

def makeHTML(sequences, annotation, VALS_VALID, styles=None, style_prefix="annot"):
    # contants
    VALS_STYLE = VALS_VALID+'?'
    
    #if VALS_STYLE==None:
    #    VALS_STYLE = VALS_VALID+'X'
    
    if None==styles:
        #('background-color:#CCCCCC;','background-color:#FFFF99;','background-color:#FF9999;'),
        #('color:#000000;','color:#00FF00;','color:#0000FF;')
        styles = []
        for i in xrange(len(annotation)):
            stl = []
            for j in xrange(len(VALS_STYLE)):
                stl.append('<span class="%s%d%s">' % (style_prefix, i, VALS_STYLE[j]))
            styles.append(stl)
    
    assert len(annotation)<=len(styles)
    for s in styles:
        assert len(s)==4
    
    allseq   = sequences+annotation
    
    output=[]
    for i in xrange(len(allseq)):
        output.append(array('c'))
    
    prevstyle=''
    dostyle=True
    for i in xrange(len(allseq[0])):
        style=''
        stylecount=0
        iann=0
        while iann<len(annotation):
            if annotation[iann][i] in VALS_VALID:
                ind = VALS_VALID.index(annotation[iann][i])
            else:
                ind = -1
            #print ind
            style+=styles[iann][ind]
            iann+=1
        
        iout=0
        dostyle = style!=prevstyle
        if dostyle:
            prevstyle=style
            prevstylecount=stylecount
        
        while iout < len(allseq):
            if dostyle:
                if i>0:
                    output[iout].extend('</span>'*len(annotation))
                output[iout].extend(style)
            output[iout].append(allseq[iout][i])
            iout+=1
    
    for iout in xrange(len(output)):
        if prevstyle:
            output[iout].extend('</span>'*len(annotation))
        output[iout] = output[iout].tostring()
    
    return output


def as_html(ali, annot_values, HTML_HEAD="", HTML_TAIL="", priority=None):
    seq   = []
    annot = []
    for eg in ali:
      if eg[0].getType() in ["structure", "sequence"]:
        seq.append(eg[0])
        annot.extend(eg[1:])
      else:
        annot.extend(eg)
    entries = seq + annot

    htmlresults = dict()
    
    for a in annot:
        # Only use the first one of every type of annotation,
        # except if we're told to prioritise a certain entry.
        # In this case, allow that to overwrite any existing
        # annotation.
        #
        if htmlresults.has_key(a.desc):
          if not priority or priority and a.code != priority:
            continue
          
        otherannots = annot[:]
        otherannots.remove(a)

        htmltext = ''
        
        htmltext += HTML_HEAD

        htmlseqs = makeHTML([e.seq for e in seq], [a.seq], annot_values)
        
        for i in xrange(len(seq)):
            htmltext += makeHTMLrow(seq[i].code, seq[i].desc, htmlseqs[i]) + "\n"
        
        htmltext += makeHTMLrow(a.code, a.desc, htmlseqs[-1]) + "\n"
        
        for a2 in otherannots:
            hlist = makeHTML([], [a2.seq], annot_values)
            htmltext += makeHTMLrow(a2.code, a2.desc, hlist[0]) + "\n"
        
        htmltext += HTML_TAIL
        
        htmlresults[a.desc] = htmltext
    
    return htmlresults
