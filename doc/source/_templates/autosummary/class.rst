{{ objname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
    :members:
    :undoc-members:
    :show-inheritance:
    
   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes
   .. autoattribute::
      :toctree:
   {% for item in attributes %}
      {{ name }}.{{ item }}
   {% endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods
   .. autofunction::
      :toctree:
   {% for item in methods %}
      {{ name }}.{{ item }}
   {% endfor %}
   {% endif %}
   {% endblock %}
