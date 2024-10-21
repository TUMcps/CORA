
getJSON('https://api.github.com/repos/TUMcps/CORA/commits?per_page=100',
    function(err, data) {
        if (err !== null) {
            alert('Something went wrong: ' + err);
        } else {
            // loop over all commits
            for (let i = 0; i < data.length; i++) {
                commit = data[i].commit;
                message = commit.message;
                messageText = message.split('\n');

                // format version
                version = messageText[0];
                // v2021 is different for several reasons...
                if (version == "CORA - v2021") {
                    break;
                }
                // special formatting for major releases (v20xx.0.0)
                if (version.endsWith("0.0")) {
                    stylestring = "color:var(--bs-primary)";
                } else {
                    stylestring = "";
                }

                // format date (long month name and numeric year and day)
                var date = new Date(commit.author.date);
                var options = { year: "numeric", month: "long", day: "numeric" };
                var formattedDate = "(" + date.toLocaleDateString("en-US", options) + ")";

                // format commit message
                formattedCommitMessage = "<ul>";
                for (let j = 2; j < messageText.length; j++) {
                    messageText_j = messageText[j].replace('-','');
                    if (messageText_j) { // text not empty
                        if (!messageText_j.startsWith('  ')) {
                            formattedCommitMessage += `<li style="list-style-type: disc"> ${messageText_j} </li>`;
                        } else {
                            formattedCommitMessage += messageText_j;
                        }
                    }
                }
                formattedCommitMessage += "</ul>";


                // only first list element is shown, others are collapsed
                if (i == 0){
                    isCollapsed = "true";
                    collapseShow = "show";
                    btnCollapsed = "";
                } else {
                    isCollapsed = "false";
                    collapseShow = "";
                    btnCollapsed = "collapsed";
                }

                container = document.getElementById("release-notes-list");
                container.innerHTML +=
                   `<li id="CORArelease_${i}" class="list-group-item px-0">
                    <a class="btn ${btnCollapsed}" data-bs-toggle="collapse" href="#collapse_CORArelease_${i}" role="button"
                        aria-expanded="${isCollapsed}"
                        aria-controls="collapse_CORArelease_${i}"
                        style="${stylestring}">
                        <div class="row">
                            <div class="col-11"><h5>${version} \u00A0 <small class="text-body-secondary">${formattedDate}</small></h5></div>
                            <div class="col-1"><span class="mr-3"></span></div>
                        </div>
                    </a>
                    <div class="collapse ${collapseShow}" id="collapse_CORArelease_${i}">
                    <div class="row">
                    <div class="col-12">
                        <p> ${formattedCommitMessage} </p>
                        <p> Show on <a href=${data[i].html_url} target="_blank"><i class="bi bi-github"></i> GitHub.</a></p>
                    </div>
                    </div>
                    </div>
                    </li>`;
            }
            // add horizontal line at the end
            container.innerHTML += `<li class="list-group-item px-0"></li>`;
        }
    });