---
title: Older Versions
custom_prev_url: "releases/0.16.0"
custom_prev_title: "Thingy"
---

{% for rel in site.oldreleases reversed %}

<h2> {{ rel.title }} (released {{ rel.date | date_to_string }}) </h2>

{{ rel.content }}

{% endfor %}


