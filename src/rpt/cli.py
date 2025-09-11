import typer

from .rpt import main

app = typer.Typer()
app.command()(main)

if __name__ == "__main__":
    app()
