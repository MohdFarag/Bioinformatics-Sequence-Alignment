import flet as ft
from utilities import *

def main(page: ft.Page):
    page.title = "Sequence Alignment"

    Sequence1 = ft.TextField(label="Enter Sequence1", autofocus=True)
    Sequence2 = ft.TextField(label="Enter Sequence2", autofocus=True)

    Sequences = ft.Column()
    DNA = ft.Checkbox(label="DNA", value=False)
    RNA = ft.Checkbox(label="RNA", value=False)
    PROTEIN = ft.Checkbox(label="PROTEIN", value=False)

    selected_files = ft.Text()

    def pick_files_result(e: ft.FilePickerResultEvent):
        selected_files.value = (
            ", ".join(map(lambda f: f.name, e.files)) if e.files else "Cancelled!"
        )
        selected_files.update()

    pick_files_dialog = ft.FilePicker(on_result=pick_files_result)

    page.overlay.append(pick_files_dialog)

    def btn_click(e):
        Sequences.controls.clear()
        Sequences.controls.append(ft.Text(f"Sequence 1 are :  {Sequence1.value} "))
        Sequences.controls.append(ft.Text(f"Sequence 2 are :  {Sequence1.value} "))

        result1, result2 = global_alignment(Sequence1.value, Sequence2.value)

        Sequences.controls.append(ft.Text(f"Sequence 1 are :  {result1} "))
        Sequences.controls.append(ft.Text(f"Sequence 2 are :  {result2} "))
        page.update()

    page.add(
        ft.Row(
            [
                ft.ElevatedButton(
                    "Pick files",
                    icon=ft.icons.UPLOAD_FILE,
                    on_click=lambda _: pick_files_dialog.pick_files(
                        allow_multiple=True)),
                          selected_files,
            ]
        ),
        ft.Row([DNA,RNA,PROTEIN]),
        Sequence1,
        Sequence2,
        ft.ElevatedButton("show Sequence", on_click=btn_click),
        Sequences,
    
    )
  
ft.app(target=main)

