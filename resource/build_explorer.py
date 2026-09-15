"""Inline the resource JSON into the explorer template.

Produces a single self-contained HTML file that works offline: no CDN, no server,
no build step at read time.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "explorer_template.html")
DATA = os.path.join(HERE, "pleiotropy_resource.json")
OUT = os.path.join(HERE, "pleiotropy_explorer.html")


def main():
    tpl = open(TEMPLATE, encoding="utf-8").read()
    raw = open(DATA, encoding="utf-8").read()
    json.loads(raw)                      # fail loudly if the data are malformed
    # a literal </script> inside the payload would close the block early
    safe = raw.replace("</", "<\\/")
    html = tpl.replace("/*__DATA__*/", safe)
    assert "/*__DATA__*/" not in html
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(html)
    print("wrote %s (%.0f KB)" % (OUT, os.path.getsize(OUT) / 1024))


if __name__ == "__main__":
    main()
