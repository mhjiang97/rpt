import typer

from .rpt import generate_report

app = typer.Typer()
app.command()(generate_report)

if __name__ == "__main__":
    app()
