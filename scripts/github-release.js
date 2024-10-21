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

    // check if confetti (<= 7 days)
    doConfetti = seconds <= 8 * 24 * 60 * 60;

    return [timeText,doConfetti]
}


getJSON('https://api.github.com/repos/TUMcps/CORA/commits?per_page=1',
    function(err, data) {
        if (err !== null) {
            alert('Something went wrong: ' + err);
        } else {
            commit = data[0].commit;
            message = commit.message;
            date = commit.author.date;
            [timeText, doConfetti] = timeSince(Date.parse(date));

            container = document.getElementById("release-version");
            container.innerHTML = `<object class="align-middle" data="https://img.shields.io/static/v1?label=Last update&message=${timeText} ago&color=4596FF" alt="TUMcps - CORA"></object>`;

            // confetti
            if (timeText, doConfetti) {
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
    });