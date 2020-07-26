---
layout: page
permalink: /teaching/
title: Teaching
description: Links to Teaching Materials - Lab Instructions, Codes, and Slides
---

<ul class="teaching-posts">

    {% assign last_year = "0" %}
    
    {% for course in site.teaching reversed %}
                
        {% if course.year != last_year %}
            
            <br>
            
            {% if course.year == 9999 %}
            
                <h3 class="year"> &nbsp; </h3>
                
            {% else %}
            
                <h3 class="year">{{ course.year }}</h3>
                
            {% endif %}
            
        {% endif %}
        
        <h4><a href="{{ course.url | prepend: site.baseurl | prepend: site.url }}">
           <strong>{{ course.title }}</strong>
        </a></h4>
        <h5>{{ course.description }}</h5>
                
        {% assign last_year = course.year %}
                
    {% endfor %}

</ul>