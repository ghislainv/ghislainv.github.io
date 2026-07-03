document.addEventListener("DOMContentLoaded", () => {
  const user = "ghislain.vieilledent";
  const domain = "cirad.fr";
  document.querySelectorAll(".email-placeholder").forEach(el => {
    el.innerHTML = '<a href="mailto:' + user + '@' + domain + '">'
                 + user + '@' + domain + '</a>';
  });
});
