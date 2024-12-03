document.addEventListener("DOMContentLoaded", function() {
    var codeElements = document.querySelectorAll('pre > code.nohighlight');
    codeElements.forEach(function(code) {
        code.parentElement.classList.add('output');
    });
});
