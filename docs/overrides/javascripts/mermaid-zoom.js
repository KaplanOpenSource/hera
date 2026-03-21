// Mermaid diagram expand button
// Polls for rendered SVGs inside .mermaid elements, then adds an Expand button.

(function () {
  function addButtons() {
    document.querySelectorAll("pre.mermaid svg, .mermaid svg").forEach(function (svg) {
      var container = svg.parentElement;
      if (!container || container.getAttribute("data-zoom-done")) return;
      container.setAttribute("data-zoom-done", "true");
      container.style.position = "relative";

      // Expand button
      var btn = document.createElement("button");
      btn.textContent = "Expand";
      btn.title = "View diagram fullscreen";
      btn.style.cssText = "position:absolute;top:4px;right:4px;z-index:10;" +
        "background:#4051b5;color:#fff;border:none;border-radius:4px;" +
        "padding:4px 10px;font-size:12px;cursor:pointer;opacity:0.8;";
      btn.onmouseenter = function () { btn.style.opacity = "1"; };
      btn.onmouseleave = function () { btn.style.opacity = "0.8"; };
      container.appendChild(btn);

      // Fullscreen overlay
      var overlay = document.createElement("div");
      overlay.style.cssText = "display:none;position:fixed;top:0;left:0;width:100vw;" +
        "height:100vh;z-index:9999;background:rgba(255,255,255,0.97);" +
        "align-items:center;justify-content:center;padding:2rem;" +
        "box-sizing:border-box;overflow:auto;cursor:zoom-out;flex-direction:column;";
      var hint = document.createElement("div");
      hint.textContent = "Click anywhere or press Escape to close";
      hint.style.cssText = "position:fixed;bottom:1rem;left:50%;transform:translateX(-50%);" +
        "font-size:14px;color:#666;";
      overlay.appendChild(hint);
      document.body.appendChild(overlay);

      // Open
      btn.addEventListener("click", function (e) {
        e.preventDefault();
        e.stopPropagation();
        overlay.querySelectorAll("svg").forEach(function (s) { s.remove(); });
        var clone = svg.cloneNode(true);
        clone.style.cssText = "max-width:95vw;max-height:85vh;width:auto;height:auto;";
        clone.removeAttribute("width");
        clone.removeAttribute("height");
        overlay.insertBefore(clone, hint);
        overlay.style.display = "flex";
        document.body.style.overflow = "hidden";
      });

      // Close
      overlay.addEventListener("click", function () {
        overlay.style.display = "none";
        document.body.style.overflow = "";
      });
    });
  }

  // Poll for Mermaid SVGs (rendered asynchronously)
  var tries = 0;
  var timer = setInterval(function () {
    addButtons();
    if (++tries > 30) clearInterval(timer);
  }, 500);

  // Escape to close all overlays
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
