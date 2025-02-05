function copyURL() {
    // Define the base domain
    const baseDomain = "https://cora.in.tum.de";

    // Get the current path, ensuring "/index.html" is removed
    let path = window.location.pathname;
    path = path.replace("/index.html", "/");
    path = path.replace("/cora_website/", "/");
    path = path.replace("/CORA/", "/");

    // Construct the final URL
    const finalURL = baseDomain + path;

    // Copy to clipboard
    navigator.clipboard.writeText(finalURL).then(() => {
        let tooltip = document.getElementById("copyTooltip");
        tooltip.innerText = "Copied!";
        setTimeout(() => { tooltip.innerText = "Copy URL"; }, 1500);
    }).catch(err => console.error("Error copying URL: ", err));
}

// load navbar
const navLinks = document.querySelectorAll('.nav-item')
const menuToggle = document.getElementById('navbarToggler')
const bsCollapse = new bootstrap.Collapse(menuToggle, {toggle: false})
navLinks.forEach((l) => {
    l.addEventListener('click', () => {
        bsCollapse.toggle()
    })
})

// add copy link to headline
headline = document.getElementById('headline')
headline.innerHTML = `${headline.innerHTML} <span id="copyURL" onclick="copyURL()" data-bs-toggle="tooltip" data-bs-title="Copy URL to clipboard"><i class="bi bi-link-45deg"></i></span>`
// load tooltips
const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]')
const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl))