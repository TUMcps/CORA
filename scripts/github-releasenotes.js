
// The GitHub API is rate limited to 60 requests/hour per IP for unauthenticated
// requests, so hitting it on every page load intermittently returns 403.
// We cache the release notes in localStorage for a week to avoid this.
const NOTES_CACHE_KEY = 'cora-release-notes';
const NOTES_CACHE_TTL = 7 * 24 * 60 * 60 * 1000; // 1 week in ms

// keep only the fields the renderer needs, so the cache stays small
function trimCommits(data) {
    return data.map(function(c) {
        return {
            commit: { message: c.commit.message, author: { date: c.commit.author.date } },
            html_url: c.html_url,
        };
    });
}

function renderReleaseNotes(data) {
    container = document.getElementById("release-notes-list");

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

        container.innerHTML +=
           `<li id="CORArelease_${i}" class="list-group-item px-0">
            <a class="btn ${btnCollapsed}" data-bs-toggle="collapse" href="#collapse_CORArelease_${i}" role="button"
                aria-expanded="${isCollapsed}"
                aria-controls="collapse_CORArelease_${i}"
                style="${stylestring}">
                <div class="row">
                    <div class="col-11"><h5>${version}   <small class="text-body-secondary">${formattedDate}</small></h5></div>
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

// try to read a fresh value from the cache first
let cachedNotes = null;
try {
    cachedNotes = JSON.parse(localStorage.getItem(NOTES_CACHE_KEY));
} catch (e) {
    cachedNotes = null;
}

if (cachedNotes && Array.isArray(cachedNotes.data) && (Date.now() - cachedNotes.fetchedAt) < NOTES_CACHE_TTL) {
    // cache hit: render without hitting the API at all
    renderReleaseNotes(cachedNotes.data);
} else {
    // cache miss or stale: fetch from GitHub
    getJSON('https://api.github.com/repos/TUMcps/CORA/commits?per_page=100',
        function(err, data) {
            // silent fail (e.g. 403 rate limit): log it but leave the list hidden
            if (err !== null || !Array.isArray(data) || data.length === 0) {
                console.warn('CORA release notes: could not fetch commits from GitHub', err, data);
                return;
            }

            const trimmed = trimCommits(data);
            try {
                localStorage.setItem(NOTES_CACHE_KEY, JSON.stringify({
                    data: trimmed,
                    fetchedAt: Date.now(),
                }));
            } catch (e) {
                // localStorage unavailable (e.g. private mode); ignore
            }

            renderReleaseNotes(trimmed);
        });
}
