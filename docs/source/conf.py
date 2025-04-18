import os
import sys

sys.path.insert(0, os.path.abspath("../../"))

project = "BGE Toolkit"
author = "Jackie Goldstein"
release = "0.0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    'sphinxcontrib.typer',
]

autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring = False

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "alabaster"

autodoc_default_options = {
    'docstring': 'short',  # This makes the signature more concise
    'no-signatures': True
}

autodoc_docstring_signature = True
autodoc_signature_format = 'short'
autodoc_typehints = 'none'

html_static_path = ['_static']  # Make sure _static exists
html_css_files = ['custom.css']


def setup(app):
    # Monkeypatch Console before sphinxcontrib.typer uses it
    import rich.console

    original_console_init = rich.console.Console.__init__

    def patched_console_init(self, *args, **kwargs):
        if "width" not in kwargs:
            kwargs["width"] = 120  # or any width that works for your CLI
        kwargs['width'] = 100
        return original_console_init(self, *args, **kwargs)

    rich.console.Console.__init__ = patched_console_init

    from sphinxcontrib.typer import TyperDirective
    from docutils import nodes


    class CustomTyperDirective(TyperDirective):
        def run(self):
            app_path = '__main__'
            result = super().run()

            def replace_in_children(node):
                for child in node.children:
                    if isinstance(child, nodes.Text) and app_path in child:
                        new_text = child.replace(app_path, "bge-toolkit")
                        child.parent.replace(child, nodes.Text(new_text))
                    elif hasattr(child, 'children'):
                        replace_in_children(child)

            for node in result:
                if isinstance(node, nodes.container):
                    replace_in_children(node)

            return result

    # Register the custom directive 'custom_typer' instead of overriding 'typer'
    app.add_directive('custom-typer', CustomTyperDirective)
