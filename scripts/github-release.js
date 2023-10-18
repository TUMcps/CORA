var getJSON = function(url, callback) {
    // https://stackoverflow.com/questions/12460378/how-to-get-json-from-url-in-javascript
    var xhr = new XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.responseType = 'json';
    xhr.onload = function() {
        var status = xhr.status;
        if (status === 200) {
            callback(null, xhr.response);
        } else {
            callback(status, xhr.response);
        }
    };
    xhr.send();
};

function timeSince(date) {
    // https://stackoverflow.com/questions/3177836/how-to-format-time-since-xxx-e-g-4-minutes-ago-similar-to-stack-exchange-site

    var seconds = Math.floor((new Date() - date) / 1000);

    var interval = seconds / 31536000;

    if (interval > 1) {
        return Math.floor(interval) + " years";
    }
    interval = seconds / 2592000;
    if (interval > 1) {
        return Math.floor(interval) + " months";
    }
    interval = seconds / 86400;
    if (interval > 1) {
        return Math.floor(interval) + " days";
    }
    interval = seconds / 3600;
    if (interval > 1) {
        return Math.floor(interval) + " hours";
    }
    interval = seconds / 60;
    if (interval > 1) {
        return Math.floor(interval) + " minutes";
    }
    return Math.floor(seconds) + " seconds";
}


getJSON('https://api.github.com/repos/TUMcps/CORA/commits?per_page=1',
    function(err, data) {
        if (err !== null) {
            alert('Something went wrong: ' + err);
        } else {
            commit = data[0].commit;
            message = commit.message;
            date = commit.author.date;

            container = document.getElementById("release-version");
            container.innerHTML = `<b>${message.split('\n')[0]}</b> <span class="text-secondary">(updated ${timeSince(Date.parse(date))} ago)</span>`;
        }
    });