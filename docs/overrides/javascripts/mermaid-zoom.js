// Click-to-expand zoom for Mermaid diagrams
document.addEventListener("DOMContentLoaded", function () {
  // MkDocs Material renders mermaid asynchronously, so we observe for new diagrams
  const observer = new MutationObserver(function () {
    document.querySelectorAll(".mermaid:not([data-zoom-bound])").forEach(function (el) {
      el.setAttribute("data-zoom-bound", "true");
      el.addEventListener("click", function (e) {
        e.stopPropagation();
        if (el.classList.contains("zoomed")) {
          el.classList.remove("zoomed");
          document.body.style.overflow = "";
        } else {
          el.classList.add("zoomed");
          document.body.style.overflow = "hidden";
        }
      });
    });
  });

  observer.observe(document.body, { childList: true, subtree: true });

  // Also close on Escape key
  document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") {
      document.querySelectorAll(".mermaid.zoomed").forEach(function (el) {
        el.classList.remove("zoomed");
        document.body.style.overflow = "";
      });
    }
  });
});
