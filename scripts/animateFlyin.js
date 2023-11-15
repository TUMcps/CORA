
// Get the elements by class name
var flyin = document.querySelectorAll(".flyin");

// Create a new intersection observer
var observer = new IntersectionObserver(function(entries) {

  // Loop through the entries
  entries.forEach(function(entry) {
    // Check if the element is intersecting
    if (entry.isIntersecting) {
      // Add the visible class
      entry.target.classList.add("visible");
    }
  });
});

// Observe the elements
flyin.forEach(function(element) {
  observer.observe(element);
});