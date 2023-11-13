function addPublicationsFromLink(listId, pubs) {
  // {bibtex: "bibtex-link", type: "pub-type"}

  // call the function to add the publication to the list
  for (let i = 0; i < pubs.length; i++){
    const pub = pubs[i];
    const bibtexlink = "../../data/bib/" + pub.bibtex;
    if (pub.type === 'diss') {
      addDissFromLink("pubsUsingCORA", bibtexlink);
    } else {
      addPaperFromLink("pubsUsingCORA", bibtexlink);
    }
  }
}

function addPaperFromLink(listId, link) {
  addPublicationFromLink(listId,link,"bi-file-text");
}

function addDissFromLink(listId, link) {
  addPublicationFromLink(listId,link,"bi-journals");
}

function addPublicationFromLink(listId, link, icon) {
  // get the unordered list element by id
  var list = document.getElementById(listId);
  // create a new list item element
  var item = document.createElement("li");

  // use the fetch API to get the bibtex entry from the link
  fetch(link)
  .then(response => response.text())
  .then(bibtex => {

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

    // check if the entry has authors
    if (entry.author) {

      // add the authors separated by commas
      author = entry.author;
      author = author.replace('\\\\\"a','&auml;')
      author = author.replace('\\\\\"o','&ouml;')
      author = author.replace('\\\\\"u','&uuml;')

      text += author;
    }

    // check if the entry has a journal
    if (entry.journal) {
      // add a period after the authors if they exist
      if (text.length > 0) {
        text += ". ";
      }
      // add the journal name in italics
      text += "<i>" + entry.journal + "</i>";
    }

    // check if the entry has a volume
    if (entry.volume) {
      // add a comma after the journal if it exists
      if (text.length > 0) {
        text += ", ";
      }
      // add the volume number
      text += entry.volume;
    }

    // check if the entry has an issue
    if (entry.issue) {
      // add a parenthesis after the volume if it exists
      if (text.length > 0) {
        text += "(";
      }
      // add the issue number
      text += entry.issue;
      // close the parenthesis
      text += ")";
    }

    // check if the entry has pages
    if (entry.pages) {
      // add a colon after the issue if it exists
      if (text.length > 0) {
        text += ": ";
      }
      // add the page range
      text += entry.pages;
    }

    // check if the entry has a year
    if (entry.year) {
      // add a period after the pages if they exist
      if (text.length > 0) {
        text += ". ";
      }
      // add the year in parentheses
      text += "(" + entry.year + ")";
    }

    // check if the entry has a url
    if (entry.url) {
      // add a period after the year if it exists
      if (text.length > 0) {
        text += ". ";
      }
      // add the url as a hyperlink
      text += "<a href='" + entry.url + "'>Link</a>";
    }
    
    // set the inner HTML of the list item to the text
    item.innerHTML = text;
    item.setAttribute('class','timeline-item mb-5');
    // append the list item to the list
    list.appendChild(item);
  })
  .catch(error => {
    // handle the error
    console.error(error);
    });
  }

