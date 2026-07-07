function timeSince(date) {
    // https://stackoverflow.com/questions/3177836/how-to-format-time-since-xxx-e-g-4-minutes-ago-similar-to-stack-exchange-site

    var seconds = Math.floor((new Date() - date) / 1000);

    var interval = seconds;
    var timeUnit = Math.floor(interval);
    const dividers = [31536000, 2592000, 86400, 3600, 60, 1];
    const names = ["year", "month", "day", "hour", "minute", "second"];

    // loop over time units
    timeText = timeUnit + " seconds"; // default to seconds
    for (let counter = 0; counter < dividers.length; counter++) {
        interval = seconds / dividers[counter];
        timeUnit = Math.floor(interval);
        if (interval >= 1) {
            if (timeUnit == 1) {
                timeText =  timeUnit + " " + names[counter];
            } else {
                timeText = timeUnit + " " + names[counter] + "s";
            }
            break
        }
    }

    // check if confetti (< 1 day)
    doConfetti = seconds <= 24 * 60 * 60;

    return [timeText,doConfetti]
}


// The GitHub API is rate limited to 60 requests/hour per IP for unauthenticated
// requests, so hitting it on every page load intermittently returns 403.
// We cache the latest commit date in localStorage for a week to avoid this.
const CACHE_KEY = 'cora-last-commit';
const CACHE_TTL = 7 * 24 * 60 * 60 * 1000; // 1 week in ms

function renderRelease(date) {
    [timeText, doConfetti] = timeSince(Date.parse(date));

    container = document.getElementById("release-version");
    container.innerHTML = `<object class="align-middle" data="https://img.shields.io/static/v1?label=Last update&message=${timeText} ago&color=4596FF" alt="TUMcps - CORA"></object>`;

    // confetti
    if (doConfetti) {
        const delay = ms => new Promise(res => setTimeout(res, ms));
        const yourFunction = async () => {
            await delay(1000);
            confetti({
                particleCount: 100,
                spread: 70,
                origin: { y: 0.3 },
            });
        };
        yourFunction();
    }
}

// try to read a fresh value from the cache first
let cached = null;
try {
    cached = JSON.parse(localStorage.getItem(CACHE_KEY));
} catch (e) {
    cached = null;
}

if (cached && cached.date && (Date.now() - cached.fetchedAt) < CACHE_TTL) {
    // cache hit: render without hitting the API at all
    renderRelease(cached.date);
} else {
    // cache miss or stale: fetch from GitHub
    getJSON('https://api.github.com/repos/TUMcps/CORA/commits?per_page=1',
        function(err, data) {
            // silent fail (e.g. 403 rate limit): log it but leave the badge hidden
            if (err !== null || !Array.isArray(data) || !data[0]) {
                console.warn('CORA last-update badge: could not fetch latest commit from GitHub', err, data);
                return;
            }

            const date = data[0].commit.author.date;
            try {
                localStorage.setItem(CACHE_KEY, JSON.stringify({
                    date: date,
                    fetchedAt: Date.now(),
                }));
            } catch (e) {
                // localStorage unavailable (e.g. private mode); ignore
            }

            renderRelease(date);
        });
}