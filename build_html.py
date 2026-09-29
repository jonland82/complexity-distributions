"""Build the responsive HTML edition from the paper's LaTeX source.

Requires Pandoc and Beautiful Soup. Run make_figures.py first to refresh the SVGs.
The generated HTML uses native MathML and needs no web connection.
"""

from __future__ import annotations

import html
import re
import subprocess
from pathlib import Path

from bs4 import BeautifulSoup


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "quicksort_complexity_distribution.tex"
OUTPUT = ROOT / "quicksort_complexity_distribution.html"
SITE_OUTPUT = ROOT / "index.html"
REPOSITORY_URL = "https://github.com/jonland82/complexity-distributions"


def pandoc_fragment(source: str) -> BeautifulSoup:
    result = subprocess.run(
        ["pandoc", "-f", "latex", "-t", "html5", "--mathml"],
        input=source,
        text=True,
        encoding="utf-8",
        capture_output=True,
        check=True,
    )
    return BeautifulSoup(result.stdout, "html.parser")


def source_metadata(source: str):
    title = re.search(r"\{\\LARGE\\bfseries\s+([^}]+)\}", source).group(1)
    author = re.search(r"\{\\normalsize\s+([^}]+)\}", source).group(1)
    date = re.search(r"\{\\small\s+([^}]+)\}", source).group(1)
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", source, re.S).group(1)
    return title, author, date, abstract


def labels_from_source(source: str):
    equations = re.findall(r"\\label\{(eq:[^}]+)\}", source)
    figures = re.findall(r"\\label\{(fig:[^}]+)\}", source)
    bibliography = re.findall(r"\\bibitem\{([^}]+)\}", source)
    return (
        {key: number for number, key in enumerate(equations, 1)},
        {key: number for number, key in enumerate(figures, 1)},
        {key: number for number, key in enumerate(bibliography, 1)},
    )


def add_equation_wrappers(soup: BeautifulSoup, equation_numbers: dict[str, int]):
    for math in list(soup.select('math[display="block"]')):
        parent = math.parent
        if parent.name == "p":
            before = soup.new_tag("p")
            while parent.contents and parent.contents[0] is not math:
                before.append(parent.contents[0].extract())
            if before.get_text(" ", strip=True):
                parent.insert_before(before)
            math.extract()
            wrapper = soup.new_tag("div", attrs={"class": "equation"})
            parent.insert_before(wrapper)
            wrapper.append(math)
            if not parent.get_text(" ", strip=True):
                parent.decompose()
        else:
            wrapper = soup.new_tag("div", attrs={"class": "equation"})
            math.wrap(wrapper)

        annotation = math.find("annotation", attrs={"encoding": "application/x-tex"})
        source = annotation.get_text() if annotation else ""
        labels = re.findall(r"\\label\{(eq:[^}]+)\}", source)
        if not labels:
            wrapper["class"].append("unnumbered")
            continue
        if len(labels) == 1:
            wrapper["id"] = labels[0]
        else:
            wrapper["id"] = labels[0] + "-group"
            for label in labels:
                marker = soup.new_tag("span", id=label, attrs={"class": "anchor-marker"})
                wrapper.insert(0, marker)
        numbers = soup.new_tag("span", attrs={"class": "equation-number"})
        for index, label in enumerate(labels):
            if index:
                numbers.append(soup.new_tag("br"))
            link = soup.new_tag("a", href="#" + label)
            link.string = f"({equation_numbers[label]})"
            numbers.append(link)
        wrapper.append(numbers)


def restore_figures(soup: BeautifulSoup, figure_numbers: dict[str, int]):
    alt = {
        "fig:fits": "Histograms of Quicksort comparison counts at n equals 64 and 256 with normal and lognormal fit curves",
        "fig:limit": "Quicksort and lognormal cumulative curves at n equals 256 and 4096, then their distinct limiting curves, with CDF gaps below",
    }
    for figure in soup.find_all("figure"):
        embed = figure.find("embed")
        if not embed:
            continue
        identifier = embed.get("id")
        image = soup.new_tag("img")
        image["src"] = Path(embed["src"]).with_suffix(".svg").name
        image["alt"] = alt.get(identifier, "Paper figure")
        image["loading"] = "lazy"
        image["decoding"] = "async"
        figure["id"] = identifier
        embed.replace_with(image)
        caption = figure.find("figcaption")
        if caption:
            caption.attrs.pop("aria-hidden", None)
            prefix = soup.new_tag("strong")
            prefix.string = f"Figure {figure_numbers[identifier]}. "
            caption.insert(0, prefix)


def restore_citations(soup: BeautifulSoup, source: str, bibliography_numbers: dict[str, int]):
    calls = re.findall(r"\\cite(?:\[([^]]+)\])?\{([^}]+)\}", source)
    spans = soup.select("span.citation[data-cites]")
    if len(calls) != len(spans):
        raise ValueError(f"Citation count differs: source={len(calls)}, HTML={len(spans)}")
    for span, (locator, keys_text) in zip(spans, calls):
        span.clear()
        span.append("[")
        for index, key in enumerate(keys_text.split(",")):
            key = key.strip()
            if index:
                span.append(", ")
            link = soup.new_tag("a", href=f"#ref-{key}")
            link.string = str(bibliography_numbers[key])
            span.append(link)
        if locator:
            span.append(", " + locator.replace("~", " ").replace("--", "–"))
        span.append("]")


def restore_references(soup: BeautifulSoup, equation_numbers, figure_numbers):
    for link in soup.select("a[data-reference-type]"):
        target = link.get("data-reference", "")
        if target.startswith("eq:"):
            link.string = f"({equation_numbers[target]})"
        elif target.startswith("fig:"):
            link.string = str(figure_numbers[target])
        link.attrs.pop("data-reference-type", None)
        link.attrs.pop("data-reference", None)


def restore_bibliography(soup: BeautifulSoup, source: str):
    old = soup.select_one(".thebibliography")
    if old:
        old.decompose()
    bibliography = re.search(
        r"\\begin\{thebibliography\}\{[^}]+\}.*?(.*?)\\end\{thebibliography\}",
        source,
        re.S,
    ).group(1)
    items = re.findall(r"\\bibitem\{([^}]+)\}(.*?)(?=\\bibitem\{|\Z)", bibliography, re.S)
    heading = soup.new_tag("h1", id="references")
    heading.string = "References"
    soup.append(heading)
    listing = soup.new_tag("ol", attrs={"class": "references"})
    for key, item in items:
        rendered = pandoc_fragment(item.strip())
        entry = soup.new_tag("li", id=f"ref-{key}")
        for child in list(rendered.contents):
            entry.append(child.extract())
        listing.append(entry)
    soup.append(listing)


def number_headings(soup: BeautifulSoup):
    sections = []
    for index, heading in enumerate(soup.find_all("h1")):
        if heading.get("id") == "references":
            heading["data-title"] = "References"
            sections.append(("references", "References", ""))
            continue
        number = str(index + 1) if index < 5 else f"Appendix {chr(ord('A') + index - 5)}"
        prefix = soup.new_tag("span", attrs={"class": "section-number", "aria-hidden": "true"})
        prefix.string = number
        heading.insert(0, prefix)
        heading["data-title"] = heading.get_text(" ", strip=True).removeprefix(number).strip()
        sections.append((heading["id"], heading["data-title"], number))
    return sections


def build():
    source = SOURCE.read_text(encoding="utf-8")
    title, author, date, abstract = source_metadata(source)
    equation_numbers, figure_numbers, bibliography_numbers = labels_from_source(source)
    soup = pandoc_fragment(source)
    if soup.select_one("div.center"):
        soup.select_one("div.center").decompose()
    keywords = soup.find("p", string=lambda text: text and text.startswith("Keywords:"))
    if not keywords:
        keywords = next((p for p in soup.find_all("p") if p.get_text().startswith("Keywords:")), None)
    keyword_html = str(keywords.extract()) if keywords else ""
    restore_figures(soup, figure_numbers)
    restore_citations(soup, source, bibliography_numbers)
    restore_references(soup, equation_numbers, figure_numbers)
    add_equation_wrappers(soup, equation_numbers)
    restore_bibliography(soup, source)
    sections = number_headings(soup)

    abstract_html = "".join(str(node) for node in pandoc_fragment(abstract).contents)
    toc = "\n".join(
        f'<a href="#{html.escape(identifier)}" data-section="{html.escape(identifier)}">'
        f'<span>{html.escape(number)}</span> {html.escape(name)}</a>'
        for identifier, name, number in sections
    )
    content = str(soup)
    page = TEMPLATE.replace("__TITLE__", html.escape(title))
    page = page.replace("__AUTHOR__", html.escape(author))
    page = page.replace("__DATE__", html.escape(date))
    page = page.replace("__ABSTRACT__", abstract_html)
    page = page.replace("__KEYWORDS__", keyword_html)
    page = page.replace("__TOC__", toc)
    page = page.replace("__PAPER__", content)
    page = page.replace("__REPOSITORY_URL__", REPOSITORY_URL)
    OUTPUT.write_text(page, encoding="utf-8")
    SITE_OUTPUT.write_text(page, encoding="utf-8")
    print(f"Wrote {OUTPUT.name} and {SITE_OUTPUT.name}: {len(sections)} sections, {len(equation_numbers)} numbered equations")


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light">
<meta name="description" content="Quicksort comparison counts, lognormal fits, and their distinct limiting distribution.">
<title>__TITLE__ — Jonathan R. Landers</title>
<style>
  :root { color-scheme: light; --ink: #242622; --muted: #5b625e; --line: #c7cbc7; --link: #205b71; --paper: #fff; --ground: #f3f3ef; }
  * { box-sizing: border-box; }
  html { scroll-behavior: auto; }
  body { margin: 0; color: var(--ink); background: var(--ground); font-family: Georgia, 'Times New Roman', serif; }
  a { color: var(--link); text-underline-offset: .15em; }
  a:hover { color: #123f51; }
  a:focus-visible, button:focus-visible { outline: 3px solid #d36d31; outline-offset: 3px; }
  .skip { position: absolute; left: -9999px; top: 1rem; z-index: 20; background: white; padding: .7rem 1rem; }
  .skip:focus { left: 1rem; }
  .layout { display: grid; grid-template-columns: minmax(180px, 235px) minmax(0, 820px); gap: clamp(1.5rem, 4vw, 4.2rem); max-width: 1220px; margin: 0 auto; padding: 1.5rem 4.5rem 5rem; }
  .contents { align-self: start; position: sticky; top: 1.4rem; max-height: calc(100vh - 2.8rem); overflow-y: auto; padding: .3rem 0; font-family: Arial, Helvetica, sans-serif; }
  .contents .toc-heading { color: #343a3b; font-size: .76rem; font-weight: 700; letter-spacing: .11em; text-transform: uppercase; padding: 0 .65rem .75rem; }
  .toc-links { display: flex; flex-direction: column; gap: .15rem; }
  .toc-links a { display: grid; grid-template-columns: 5.3rem 1fr; gap: .35rem; padding: .56rem .65rem; border-left: 2px solid transparent; color: #4c5657; text-decoration: none; font-size: .82rem; line-height: 1.33; border-radius: 0 4px 4px 0; }
  .toc-links a span { font-weight: 700; white-space: nowrap; }
  .toc-links a:hover, .toc-links a.active { background: #e8edeb; color: #153c48; border-left-color: #2b6677; }
  .paper { min-width: 0; background: var(--paper); box-shadow: 0 2px 24px #2b302216; padding: clamp(1.3rem, 4vw, 3.4rem) clamp(1.2rem, 5vw, 4.2rem) 5rem; }
  .running-head { display: flex; justify-content: space-between; flex-wrap: wrap; gap: .5rem 1rem; border-bottom: 1px solid #555; padding-bottom: .45rem; color: #414541; font: .68rem Arial, Helvetica, sans-serif; letter-spacing: .08em; text-transform: uppercase; }
  .running-head-links { display: flex; flex-wrap: wrap; gap: .4rem .9rem; }
  .running-head a { text-transform: none; letter-spacing: 0; white-space: nowrap; }
  .hero { text-align: center; padding: 2rem 0 .7rem; }
  .hero h1 { border: 0; margin: 0 auto 1rem; max-width: 680px; font-size: clamp(1.8rem, 3vw, 2.35rem); line-height: 1.14; letter-spacing: -.022em; }
  .hero .author { margin: 0; font-size: 1.05rem; }
  .hero .date { margin: .4rem 0 0; font-size: .85rem; color: var(--muted); }
  .abstract { max-width: 700px; margin: 1.3rem auto .8rem; font-size: .93rem; line-height: 1.5; }
  .abstract h2 { margin: 0 0 .55rem; text-align: center; font-size: .92rem; }
  .abstract p { margin: 0; }
  .keywords { font-size: .86rem; margin: 1.05rem 0 1.8rem; }
  .paper-body { font-size: 1.01rem; line-height: 1.59; }
  .paper-body h1 { scroll-margin-top: 5.5rem; font-size: clamp(1.35rem, 2vw, 1.58rem); line-height: 1.25; margin: 2.45rem 0 .9rem; padding: 0 0 .35rem; border-bottom: 1px solid #5d625c; }
  .section-number { display: inline-block; margin-right: .75rem; }
  .paper-body p { margin: .7rem 0; }
  .paper-body figure { margin: 1.7rem 0 2rem; scroll-margin-top: 5rem; }
  .paper-body figure img { display: block; width: 100%; height: auto; }
  .paper-body figcaption { margin-top: .7rem; font-size: .84rem; line-height: 1.45; color: #484d4a; }
  .paper-body table { border-collapse: collapse; margin: 1rem auto; font-size: .92rem; }
  .paper-body td, .paper-body th { padding: .25rem .65rem; text-align: center; }
  .paper-body tr:first-child { border-top: 1px solid #484d4a; }
  .paper-body tr:last-child { border-bottom: 1px solid #484d4a; }
  .paper-body .center { text-align: center; overflow-x: auto; }
  .paper-body math[display="inline"] { font-size: 1.025em; }
  .equation { position: relative; display: grid; grid-template-columns: minmax(0, 1fr) max-content minmax(0, 1fr); align-items: center; width: 100%; overflow-x: auto; padding: .65rem .1rem; margin: .2rem 0 .55rem; scroll-margin-top: 5rem; }
  .equation math { grid-column: 2; font-size: 1.05em; }
  .equation-number { grid-column: 3; justify-self: end; min-width: 2.3rem; text-align: right; font-size: .83rem; line-height: 1.8; }
  .equation-number a { color: #4c5250; text-decoration: none; }
  .anchor-marker { position: absolute; top: 0; width: 0; height: 0; scroll-margin-top: 5rem; }
  .paper-body .lemma, .paper-body .theorem { padding: .15rem .9rem; margin: 1.3rem 0; border-left: 3px solid #4e7880; background: #f6f8f7; scroll-margin-top: 5rem; }
  .paper-body .proof { margin: .85rem 0 1.4rem; }
  .paper-body .references { padding-left: 1.5rem; font-size: .88rem; overflow-wrap: anywhere; }
  .paper-body .references li { padding-left: .3rem; margin: .55rem 0; scroll-margin-top: 5rem; }
  .paper-body .references p { margin: 0; }
  .paper-footer { border-top: 1px solid var(--line); margin-top: 2.6rem; padding-top: .8rem; color: var(--muted); font: .8rem Arial, Helvetica, sans-serif; }
  .paper-body .citation { white-space: nowrap; font-size: .9em; }
  .paper-body code { font-size: .85em; background: #f2f3f1; padding: .05em .24em; }
  .edge-nav { position: fixed; z-index: 5; top: 50%; transform: translateY(-50%); display: flex; align-items: center; justify-content: center; width: 3rem; height: 4.5rem; background: #fffef9e8; border: 1px solid #ccd1cc; color: #285a66; text-decoration: none; font: 1.8rem Arial, sans-serif; box-shadow: 0 3px 12px #00000012; }
  .edge-nav:hover { background: white; }
  .edge-nav.previous { left: .5rem; border-radius: 5px; }
  .edge-nav.next { right: .5rem; border-radius: 5px; }
  .edge-nav.disabled { pointer-events: none; opacity: .28; }
  .mobile-toc, .mobile-arrows { display: none; }
  @media (max-width: 1050px) { .layout { padding-left: 3.8rem; padding-right: 3.8rem; gap: 1.1rem; grid-template-columns: 170px minmax(0, 1fr); } .toc-links a { grid-template-columns: 1fr; gap: .05rem; } }
  @media (max-width: 760px) {
    .layout { display: block; padding: 0 0 4.2rem; }
    .contents, .edge-nav { display: none; }
    .paper { box-shadow: none; padding: 1.3rem clamp(1rem, 4.4vw, 2rem) 4.5rem; }
    .mobile-toc { position: sticky; top: 0; z-index: 8; display: flex; gap: .35rem; overflow-x: auto; white-space: nowrap; padding: .5rem .7rem; background: #f5f5f1f5; border-bottom: 1px solid var(--line); font: .76rem Arial, sans-serif; scrollbar-width: thin; }
    .mobile-toc::-webkit-scrollbar { height: 3px; }
    .mobile-toc a { color: #36545b; text-decoration: none; padding: .42rem .7rem; border-radius: 999px; background: #e8edeb; }
    .mobile-toc a.active { background: #295d6c; color: white; }
    .mobile-arrows { position: fixed; bottom: 0; left: 0; right: 0; z-index: 8; display: flex; justify-content: space-between; gap: 1rem; padding: .55rem max(.8rem, env(safe-area-inset-left)) calc(.55rem + env(safe-area-inset-bottom)) max(.8rem, env(safe-area-inset-right)); background: #fffdf9f5; border-top: 1px solid var(--line); font: .82rem Arial, sans-serif; }
    .mobile-arrows a { text-decoration: none; padding: .4rem .6rem; }
    .mobile-arrows a.disabled { opacity: .3; pointer-events: none; }
    .running-head { font-size: .57rem; }
    .hero { padding-top: 1.65rem; }
    .paper-body { font-size: .97rem; }
    .paper-body h1 { scroll-margin-top: 4.8rem; }
    .equation math { font-size: .97em; }
  }
  @media print { body { background: white; } .layout { display: block; max-width: none; padding: 0; } .paper { box-shadow: none; padding: 0; } .contents, .mobile-toc, .edge-nav, .mobile-arrows, .running-head a { display: none; } .paper-body h1, figure, .equation { break-inside: avoid; } a { color: inherit; } }
</style>
</head>
<body>
<a class="skip" href="#paper-content">Skip to paper</a>
<nav class="mobile-toc" aria-label="Paper sections">__TOC__</nav>
<a class="edge-nav previous disabled" id="edge-prev" href="#top" aria-label="Previous section" title="Previous section">‹</a>
<a class="edge-nav next" id="edge-next" href="#top" aria-label="Next section" title="Next section">›</a>
<div class="layout">
  <aside class="contents" aria-label="Contents"><div class="toc-heading">Contents</div><nav class="toc-links">__TOC__</nav></aside>
  <article class="paper" id="paper-content">
    <div class="running-head"><span>Landers / The Shape of Quicksort Cost</span><span class="running-head-links"><a href="quicksort_complexity_distribution.pdf">Download PDF</a><a href="__REPOSITORY_URL__">GitHub repository</a></span></div>
    <header class="hero" id="top"><h1>__TITLE__</h1><p class="author">__AUTHOR__</p><p class="date">__DATE__</p></header>
    <section class="abstract" aria-labelledby="abstract-title"><h2 id="abstract-title">Abstract</h2>__ABSTRACT__</section>
    <div class="keywords">__KEYWORDS__</div>
    <div class="paper-body">__PAPER__</div>
    <footer class="paper-footer">Source, figures, and build instructions: <a href="__REPOSITORY_URL__">GitHub repository</a>.</footer>
  </article>
</div>
<nav class="mobile-arrows" aria-label="Previous and next sections"><a id="mobile-prev" class="disabled" href="#top">← Previous</a><a id="mobile-next" href="#top">Next →</a></nav>
<script>
(() => {
  const headings = [...document.querySelectorAll('.paper-body > h1')];
  const links = [...document.querySelectorAll('[data-section]')];
  const previous = [document.getElementById('edge-prev'), document.getElementById('mobile-prev')];
  const next = [document.getElementById('edge-next'), document.getElementById('mobile-next')];
  function update() {
    const threshold = window.innerWidth <= 760 ? 110 : 85;
    let index = -1;
    for (let i = 0; i < headings.length; i++) {
      if (headings[i].getBoundingClientRect().top <= threshold) index = i;
    }
    const current = index < 0 ? null : headings[index].id;
    for (const link of links) link.classList.toggle('active', link.dataset.section === current);
    for (const link of previous) {
      const target = index <= 0 ? '#top' : '#' + headings[index - 1].id;
      link.href = target;
      link.classList.toggle('disabled', index < 0);
      link.title = index < 0 ? 'Start of paper' : 'Previous: ' + headings[index - 1].dataset.title;
    }
    for (const link of next) {
      const target = index + 1 < headings.length ? '#' + headings[index + 1].id : '#references';
      link.href = target;
      link.classList.toggle('disabled', index >= headings.length - 1);
      link.title = index + 1 < headings.length ? 'Next: ' + headings[index + 1].dataset.title : 'End of paper';
    }
  }
  document.addEventListener('scroll', update, {passive: true});
  window.addEventListener('resize', update, {passive: true});
  window.addEventListener('hashchange', update);
  update();
})();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    build()
