---
title: Complete history (most recent first)
---

{% for rel in site.releases reversed %}
  {% if rel.number == 'changes' %}
  <h2> {{ rel.title }} </h2>
  {% else %}
  <h2> {{ rel.title }} (released {{ rel.date | date_to_string }}) </h2>
  {% endif %}
  {{ rel.content }} 
{% endfor %}


