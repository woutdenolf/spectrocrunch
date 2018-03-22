{{ fullname | escape | underline }}

Description
-----------

.. automodule:: {{ fullname | escape }}
    :members:
    :undoc-members:
    :show-inheritance:

{% if classes %}
Classes
-------
.. autosummary:
    :toctree: _autosummary
    {% for class in classes %}
        {{ class }}
    {% endfor %}
{% endif %}

{% if functions %}
Functions
---------
.. autosummary:
    :toctree: _autosummary
    {% for function in functions %}
        {{ function }}
    {% endfor %}
{% endif %}
