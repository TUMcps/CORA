async function addPublicationsFromLink(listId, pubs) {
  // {bibtex: "bibtex-link", type: "pub-type", pdf: "pdf-link"}

  // call the function to add the publication to the list
  for (let i = 0; i < pubs.length; i++){
    const pub = pubs[i];
    const bibtexlink = "../../data/bib/" + pub.bibtex;
    const pdflink = pub.pdf;

    // add publication; await required to maintain order
    if (pub.type === 'diss') {
      await addDissFromLink(listId, bibtexlink, pdflink);
    } else {
      await addPaperFromLink(listId, bibtexlink, pdflink);
    }
  }
}

async function addPaperFromLink(listId, link, pdflink) {
  await addPublicationFromLink(listId,link,pdflink,"bi-file-text");
}

async function addDissFromLink(listId, link, pdflink) {
  await addPublicationFromLink(listId,link,pdflink,"bi-journals");
}

async function addPublicationFromLink(listId, link, pdflink, icon) {
  //async required due to await during bibtex fetch

  // get the unordered list element by id
  var list = document.getElementById(listId);
  // create a new list item element
  var item = document.createElement("li");

  // use the fetch API to get the bibtex entry from the link
  const response = await fetch(link);
  // get the bibtex text from the response
  const bibtex = await response.text();

  // use the bibtex-parser-js library to parse the bibtex entry
  //var parser = require('bibtex-parser-js');
  //var parsed = parser.toJSON(bibtex);
  const parsed = BibtexParser.parseToJSON(bibtex);
  // get the first entry from the parsed object
  var entry = parsed[0];
  // create a string to store the publication details
  // start with icon
  var text = "<span class=\"timeline-icon\"><i class=\"bi " + icon + "\"></i></span>";

  // check if the entry has a title
  if (entry.title) {
    // add the title in bold
    text += "<h5>" + entry.title + "</h5>";
  }


  // add date
  date = "";
  if (entry.month) {
      date += entry.month;
  }
  if (entry.year) {
    date += " " + entry.year;
  }
  text += "<p class=\"text-muted mb-2 fw-bold\"><i class=\"bi bi-calendar-event\"></i> " + date + "</p>";

  // start paragraph element to style the text
  text += "<p class=\"text-muted mb-2\"> ";
  // style journal: Author(s). Journal Volume (Issue). Pages. (Year)
  // style journal: Author(s). Conference. Pages. (Year)
  // style dissertation: Author. Dissertation. School. Pages. (Year)

  // check if the entry has authors
  if (entry.author) {

    // add the authors separated by commas
    author = entry.author;
    author = author.replaceAll('\\\\\"a','&auml;')
    author = author.replaceAll('\\\\\"o','&ouml;')
    author = author.replaceAll('\\\\\"u','&uuml;')

    // handle commas and 'and's in authors
    ands = author.match(/and/g);
    if (ands === null){
      // don't change anything
    } else if (ands.length == 1) {
      // don't change anything
    } else if (ands.length > 1) {
      // replace all but last one by ", " and add ", " before last "and"
      author = author.replaceAll(' and', ', ');
      idxLastComma = author.lastIndexOf(',');
      author = author.slice(0,idxLastComma+1) + " and" + author.slice(idxLastComma+1);
    }

    // add author to text
    text += author + ". ";
  }

  // journal?
  if (entry.journal) {
    // add the journal name in italics
    text += "<i>" + entry.journal;

    // check if the entry has a volume
    if (entry.volume) {
      // add the volume number
      text += " " + entry.volume;
    }

    // check if the entry has an issue
    if (entry.issue) {
      // add the issue number
      text += "(" + entry.issue + ")";
    } else if (entry.number) {
      // add the issue number
      text += "(" + entry.number + ")";
    }

    // end italics
    text += "</i>";
  } else if (entry.booktitle) {
    // add the proceedings name in italics
    text += "<i>" + entry.booktitle + "</i>";
  } else if (entry.school) {
    // add 'Dissertation' and the school
    text += "Dissertation, " + entry.school + ".";
  }

  // check if the entry has pages
  if (entry.pages) {
    // add the page range
    text += ". pp. " + entry.pages + ".";
  } else if (entry.articleno) {
    // sometimes, we have article number + number of pages
    text += ". Article No. " + entry.articleno;
    if (entry.numpages) {
      text += ", " + entry.numpages + " pages.";
    }
  }

  // close paragraph with author, title, etc.
  if (!text.trimEnd().endsWith('.')) {
    text += '.'
  }
  text += "</p>";

  // new div for clickable data: .bib, .pdf, URL
  text += "<div class=\"btn-group\" role=\"group\"> ";
  buttonclass = "class=\"btn btn-outline-primary btn-sm\"";

  // each entry has a .bib file
  text += "<a href='" + link + "' " + buttonclass + ">.bib</a>";

  // check if PDF exists
  if (pdflink) {
    text += "<a target=\"_blank\" href='" + pdflink + "' " + buttonclass + ">PDF</a>";
  }

  // check if the entry has a url
  if (entry.url) {

    // add the url as a hyperlink
    text += "<a target=\"_blank\" href='" + entry.url + "' " + buttonclass + ">URL</a>";
  }

  // close paragraph
  text += "</div>";

  // pill badges for keywords
  if (entry.keywords) {
    // badge class (ms-x is vertical separation)
    badgeclass = "\"badge border border-secondary text-secondary ms-2\"";

    // separate keywords into single badges
    keywords = entry.keywords.split(",");
    for (let i = 0; i < keywords.length; i++){
      text += "<div class=" + badgeclass + ">" + keywords[i].trim() + "</div>";
    }
  }

  // set the inner HTML of the list item to the text
  item.innerHTML = text;
  item.setAttribute('class','timeline-item mb-5');
  // append the list item to the list
  list.appendChild(item);
}

