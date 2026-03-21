// Mermaid diagram expand button
// Polls DOM for SVGs rendered by Mermaid and adds an Expand button.

(function () {
  function addButtons() {
    // Find all SVGs on the page that look like Mermaid diagrams
    document.querySelectorAll("svg").forEach(function (svg) {
      // Skip if already processed
      if (svg.getAttribute("data-zoom-done")) return;
      // Only target Mermaid diagrams (they have aria-roledescription or class containing mermaid-related ids)
      var isMermaid = svg.getAttribute("aria-roledescription") === "classDiagram" ||
                      svg.getAttribute("aria-roledescription") === "flowchart-v2" ||
                      svg.getAttribute("aria-roledescription") === "flowchart" ||
                      svg.getAttribute("aria-roledescription") === "sequence" ||
                      svg.getAttribute("aria-roledescription") === "stateDiagram" ||
                      svg.getAttribute("aria-roledescription") ||
                      (svg.id && svg.id.startsWith("mermaid-")) ||
                      (svg.parentElement && svg.parentElement.classList.contains("mermaid"));
      if (!isMermaid) return;

      svg.setAttribute("data-zoom-done", "true");

      var container = svg.parentElement;
      if (!container) return;
      container.style.position = "relative";

      // Expand button
      var btn = document.createElement("button");
      btn.textContent = "Expand";
      btn.title = "View diagram fullscreen";
      Object.assign(btn.style, {
        position: "absolute", top: "4px", right: "4px", zIndex: "10",
        background: "#4051b5", color: "#fff", border: "none", borderRadius: "4px",
        padding: "4px 10px", fontSize: "12px", cursor: "pointer", opacity: "0.8"
      });
      btn.onmouseenter = function () { btn.style.opacity = "1"; };
      btn.onmouseleave = function () { btn.style.opacity = "0.8"; };
      container.appendChild(btn);

      // Fullscreen overlay
      var overlay = document.createElement("div");
      Object.assign(overlay.style, {
        display: "none", position: "fixed", top: "0", left: "0",
        width: "100vw", height: "100vh", zIndex: "9999",
        background: "rgba(255,255,255,0.97)", alignItems: "center",
        justifyContent: "center", padding: "2rem", boxSizing: "border-box",
        overflow: "auto", cursor: "zoom-out", flexDirection: "column"
      });
      var hint = document.createElement("div");
      hint.textContent = "Click anywhere or press Escape to close";
      Object.assign(hint.style, {
        position: "fixed", bottom: "1rem", left: "50%",
        transform: "translateX(-50%)", fontSize: "14px", color: "#666"
      });
      overlay.appendChild(hint);
      document.body.appendChild(overlay);

      btn.addEventListener("click", function (e) {
        e.preventDefault();
        e.stopPropagation();
        overlay.querySelectorAll("svg").forEach(function (s) { s.remove(); });
        var clone = svg.cloneNode(true);
        clone.removeAttribute("width");
        clone.removeAttribute("height");
        clone.removeAttribute("style");
        Object.assign(clone.style, {
          maxWidth: "95vw", maxHeight: "85vh", width: "auto", height: "auto"
        });
        overlay.insertBefore(clone, hint);
        overlay.style.display = "flex";
        document.body.style.overflow = "hidden";
      });

      overlay.addEventListener("click", function () {
        overlay.style.display = "none";
        document.body.style.overflow = "";
      });
    });
  }

  // Poll every 500ms for up to 20 seconds (Mermaid renders async)
  var tries = 0;
  var timer = setInterval(function () {
    addButtons();
    if (++tries > 40) clearInterval(timer);
  }, 500);

  // Escape to close
  document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") {
      document.querySelectorAll("div").forEach(function (el) {
        if (el.style.zIndex === "9999" && el.style.display === "flex") {
          el.style.display = "none";
          document.body.style.overflow = "";
        }
      });
    }
  });
})();
