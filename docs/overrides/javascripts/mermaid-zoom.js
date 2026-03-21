// Click-to-expand zoom for Mermaid diagrams
// Adds an "Expand" button to each diagram and a fullscreen overlay on click.

(function () {
  function addZoomToMermaid() {
    // MkDocs Material renders mermaid as <pre class="mermaid"> then replaces with SVG
    // Target both <pre class="mermaid"> and any element with class "mermaid" containing an SVG
    document.querySelectorAll("pre.mermaid, .mermaid").forEach(function (el) {
      if (el.getAttribute("data-zoom-bound")) return;
      if (!el.querySelector("svg") && el.tagName !== "PRE") return;
      el.setAttribute("data-zoom-bound", "true");

      // Wrap in a container
      var wrapper = document.createElement("div");
      wrapper.className = "mermaid-zoom-wrapper";
      el.parentNode.insertBefore(wrapper, el);
      wrapper.appendChild(el);

      // Add expand button
      var btn = document.createElement("button");
      btn.className = "mermaid-zoom-btn";
      btn.textContent = "Expand";
      btn.title = "Click to view diagram fullscreen";
      wrapper.appendChild(btn);

      // Create overlay
      var overlay = document.createElement("div");
      overlay.className = "mermaid-zoom-overlay";
      var hint = document.createElement("div");
      hint.className = "mermaid-zoom-close-hint";
      hint.textContent = "Click anywhere or press Escape to close";
      overlay.appendChild(hint);
      document.body.appendChild(overlay);

      // Open overlay
      btn.addEventListener("click", function (e) {
        e.stopPropagation();
        var svg = el.querySelector("svg");
        if (svg) {
          var clone = svg.cloneNode(true);
          // Remove any existing SVG in overlay (keep hint)
          overlay.querySelectorAll("svg").forEach(function (s) { s.remove(); });
          overlay.insertBefore(clone, hint);
        }
        overlay.classList.add("active");
        document.body.style.overflow = "hidden";
      });

      // Close overlay on click
      overlay.addEventListener("click", function () {
        overlay.classList.remove("active");
        document.body.style.overflow = "";
      });
    });
  }

  // Run on load and observe for async Mermaid rendering
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", function () {
      setTimeout(addZoomToMermaid, 1000);
    });
  } else {
    setTimeout(addZoomToMermaid, 1000);
  }

  // Re-run on navigation (MkDocs instant loading)
  var observer = new MutationObserver(function () {
    setTimeout(addZoomToMermaid, 500);
  });
  observer.observe(document.body, { childList: true, subtree: true });

  // Close on Escape
  document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") {
      document.querySelectorAll(".mermaid-zoom-overlay.active").forEach(function (ov) {
        ov.classList.remove("active");
        document.body.style.overflow = "";
      });
    }
  });
})();
