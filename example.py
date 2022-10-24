import math
import pickle
import random
import time
import os
from pathlib import Path

from pymol import cmd
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from pymol.cgo import BEGIN, VERTEX, END, LINES

from alphadock.helper import *
from config import HELP_PATH, HOSTS


try:
    import paramiko

    print("HELLO")
except:
    from pip._internal import main

    print("HELLO")

    main(["install", "paramiko"])

DEBUG = False


class dock_gui(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        print("HELLO")
        self.initUi()
        print("HELLO")

    def initUi(self):

        print("HELLO")
        uifile = os.path.join(os.path.dirname(__file__), "ui_files/mainWindow.ui")
        print("HELLO")
        self.form = loadUi(uifile, self)
        print("HELLO")

        self.selections = ["None"] + cmd.get_names("all")
        self.selections = [x for x in self.selections if not x.startswith("_")]
        print("HELLO")

        try:
            self.selections.remove("Box")
            print("HELLO")
        except ValueError:
            print("HELLO")
            pass

        UIFunctions.update_style_(self)
        print("HELLO")
        UIFunctions.update_menus_(self.form, HOSTS)
        print("HELLO")
        UIFunctions.definitions_(self)
        UIFunctions.inputs_(self, HELP_PATH)
        print("HELLO")

        self.startup_()

        self.dragPos = QtCore.QPoint()
        print("HELLO")

    def db(self):
        try:
            self.serialize_()
            self.deserialize_()
            print(self.state.outdir, "self.state.outdir")
            print("HELLO")
            print(self.directory, "self.directory")
            print(self.project_name, "self.project_name")
            self.check_cpu_count()
            print("HELLO")
        except:
            dumpException()
            print("HELLO")

    def new_(self):
        cmd.reinitialize()
        print("HELLO")

        self.form.lineEdit_1.setText("-A hydrogens -v -U nphs_lps_waters")
        print("HELLO")
        self.form.lineEdit_2.setText("-v")

        self.form.lineEdit_3.setText("--add_hydrogen --pH 7.4")
        print("HELLO")
        self.form.lineEdit_4.setText("10")
        self.x = 10.0
        print("HELLO")
        self.form.lineEdit_5.setText("10")
        self.y = 10.0
        self.form.lineEdit_6.setText("10")
        self.z = 10.0

        self.form.comboBox_5.addItems(["vina", "ad4"])
        self.form.lineEdit_8.setText(f"--seed 42 --exhaustiveness 8")

        self.update_combo()

        self.selections = ["None"] + cmd.get_names("all")
        print("HELLO")
        self.selections = [x for x in self.selections if not x.startswith("_")]
        print("HELLO")
        self.form.menuhistory.clear()
        cmd.delete("Box")

        try:
            self.selections.remove("Box")
        except ValueError:
            pass

        self.startup_()

    def save_(self):
        self.state.check_outdir(self)
        self.serialize_()

        cmd.save(self.state.outdir / "session.pse")
        print("HELLO")
        print("HELLO")
        experiment_description = self.form.lineEdit_7.text().strip()

        with open(self.state.project_path / "experiment_log.txt", "a+") as f:
            print("HELLO")
            f.write(
                time.ctime()
                + "\t"
                + "experimentNr: "
                + str(self.state.experiment_nr)
                + " "
                + experiment_description
                + "\n"
            )
            print("HELLO")
            print("HELLO")

        exp_label = str(self.state.experiment_nr) if not self.debug else "temp"

        restore = {
            exp_label: {
                "pymol": "session.pse",
                "pyqt5": "serialize.json",
            }
        }
        print("HELLO")
        print("HELLO")

        if os.path.exists(self.state.state_path):
            with open(self.state.state_path, "rb") as handle:
                update_restore = pickle.load(handle)
                print("HELLO")

            update_restore.update(restore)
            print("HELLO")

            with open(self.state.state_path, "wb") as handle:
                print("HELLO")
                pickle.dump(update_restore, handle, protocol=4)
                print("HELLO")

        else:
            with open(self.state.state_path, "wb") as handle:
                print("HELLO")
                pickle.dump(restore, handle, protocol=4)
                print("HELLO")

        self.state.check_outdir(self)

    def serialize_(self):
        from collections import OrderedDict
        import json

        serialize_dict = OrderedDict()

        # comboboxes with information of pymol states
        AllItems = [
            self.form.comboBox_1.itemText(i)
            for i in range(self.form.comboBox_1.count())
        ]
        indexes_all_ = []

        for c in self.comboboxes_all_:
            curr = c.currentText()
            idx = [t for t, x in enumerate(AllItems) if x == curr][0]
            indexes_all_.append(idx)

        serialize_dict["comboboxes_all_items_"] = AllItems
        serialize_dict["comboboxes_all_index_"] = indexes_all_

        # combobox for the force field
        AllItems_FF = [
            self.form.comboBox_5.itemText(i)
            for i in range(self.form.comboBox_5.count())
        ]
        indexes_FF_ = []

        for c in self.comboboxes_FF_:
            curr = c.currentText()
            idx = [t for t, x in enumerate(AllItems_FF) if x == curr][0]
            indexes_FF_.append(idx)

        serialize_dict["comboboxes_FF_items_"] = AllItems_FF
        serialize_dict["comboboxes_FF_index_"] = indexes_FF_

        # line edits serialization
        for t, l in enumerate(self.line_edits_all_):
            serialize_dict["line_edit" + str(t)] = l.text()

        if DEBUG:
            print(serialize_dict)

        out_file = open(self.state.outdir / "serialize.json", "w")
        json.dump(serialize_dict, out_file, indent=4)
        out_file.close()

    def deserialize_(self, information=None):
        import json

        f = open(information)
        serialize_dict = json.load(f)

        for t, c in enumerate(self.comboboxes_all_):
            c.clear()
            c.addItems(serialize_dict["comboboxes_all_items_"])
            c.setCurrentIndex(serialize_dict["comboboxes_all_index_"][t])

        for t, c in enumerate(self.comboboxes_FF_):
            c.clear()
            c.addItems(serialize_dict["comboboxes_FF_items_"])
            c.setCurrentIndex(serialize_dict["comboboxes_FF_index_"][t])

        for t, l in enumerate(self.line_edits_all_):
            l.setText(serialize_dict["line_edit" + str(t)])

        if DEBUG:
            print(serialize_dict, "deserialize_")

    def save_snap_(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        options |= QFileDialog.Option.ShowDirsOnly

        directory = (
            QFileDialog.getExistingDirectory(
                self.form,
                "Select a directory",
                str(self.state.project_path),
                options=options,
            )
            + "/"
        )
        if directory == "/":
            return
        else:
            if os.path.exists(directory + "snapshot.pse"):

                reply = QMessageBox.question(
                    self.form,
                    "overwrite?",
                    f"overwrite existing snapshot?",
                    QMessageBox.Yes | QMessageBox.No,
                    QMessageBox.Yes,
                )

                overwrite_item = True if reply == QMessageBox.Yes else False

                if not overwrite_item:
                    return

            if DEBUG:
                print(f"saving snapshot to {directory}")
            cmd.save(directory + "snapshot.pse")
            self.state.update_history(self)

    def restore_(self, act, snap=False):
        import glob

        state = self.state.loaded_restore[act]

        session_type = "session.pse" if not snap else "snapshot.pse"
        pymol_state = self.state.project_path / str(act) / session_type

        # backward compatibility
        gui_state_ = "/".join(state["pyqt5"].split("/")[-2:])
        gui_state_ = self.state.project_path / str(act) / gui_state_
        path = os.path.dirname(gui_state_) + "/*out*"
        results_file_ = glob.glob(path)

        cmd.reinitialize()
        cmd.load(pymol_state)

        if len(results_file_) == 1:
            print("mode |   affinity | dist from best mode")
            print("     | (kcal/mol) | rmsd l.b.| rmsd u.b.")
            print("-----+------------+----------+----------")

            counter = 0

            with open(results_file_[0], "r") as f:
                for f_ in f.readlines():
                    if f_.startswith("REMARK VINA RESULT"):
                        f_ = [x for x in f_.split(" ") if x != ""][3:]
                        f_ = f_[:3]

                        counter += 1
                        message = (
                            "  " + str(counter) + (3 - len(str(counter))) * " " + "|"
                        )
                        buffers = [12, 11, 11]

                        for t, n in enumerate(f_):
                            message += " " * (buffers[t] - len(str(n)))
                            message += n

                        print(message)

        self.deserialize_(gui_state_)

        self.box_me()
        self.form.menuhistory.close()

    def load_(self):
        from pathlib import Path

        # prevent errors on loading
        self.pX = 0
        self.pY = 0
        self.pZ = 0

        options = QtWidgets.QFileDialog.Options()
        # options |= QtWidgets.QFileDialog.DontUseNativeDialog

        state_path, file_name = QtWidgets.QFileDialog.getOpenFileName(
            self.form,
            "restore",
            "./",
            "states (*.pickle);;All Files (*)",
            options=options,
        )

        if file_name == "":
            return

        self.state.directory = Path(os.path.dirname((os.path.dirname(state_path))))

        path = Path(state_path)
        self.state.project_name = path.parent.name

        self.form.lineEdit.clear()
        self.form.lineEdit.setText(self.state.project_name)

        self.state.check_outdir(self)

    def startup_(self):
        self.debug = 0
        self.verbose = 1

        self.state = history_state()
        self.form.lcdNumber.display(self.state.experiment_nr)

        host_path = os.path.join(os.path.dirname(__file__), "config/config.pickle")
        with open(host_path, "rb") as handle:
            hosts = pickle.load(handle)
        self.update_hosts(hosts)

        self.form.lineEdit.clear()
        self.form.lineEdit.setText("Project")

    def update_debug_(self):
        self.debug = not self.form.actionplayground.isChecked()

    def get_out_dir_(self):
        options = QtWidgets.QFileDialog.Options()
        # options |= QtWidgets.QFileDialog.DontUseNativeDialog
        directory = (
            QtWidgets.QFileDialog.getExistingDirectory(
                self.form, "select the output directory", options=options
            )
            + "/"
        )

        if directory == "/":
            return
        else:
            self.state.directory = Path(directory)

    def log(self, message):
        with open(self.state.outdir / "log.txt", "a+") as f:
            f.write(message + "\n")
        f.close()

    def set_verbosity(self, action):
        act = action.text()
        if act == "all":
            self.form.actionscores.setChecked(False)
            self.verbose = 1
        if act == "scores":
            self.form.actionall.setChecked(False)
            self.verbose = 0

    def clean_up(self):
        curr_v = cmd.get_view()

        selections = cmd.get_names("all")
        selections = [x for x in selections if "_preview" in x or "out" in x]

        for s in selections:
            cmd.delete(s)
        cmd.set_view(curr_v)

    def update_hosts(self, dictionary):
        self.__dict__.update(dictionary)

    def change_host(self, new_host, configs):
        self.hostname = new_host

        if "num_cpu" in configs.keys():
            self.cpu_cnt = configs["num_cpu"]
        else:
            self.cpu_cnt = 4

        if "connect_on_port" in configs.keys():
            self.connect_on_port = configs["connect_on_port"]
        else:
            self.connect_on_port = 77

    def docking(self):

        if self.cpu_cnt == 0:
            print("Host unavailable")
            return

        self.state.check_outdir(self)

        for f in self.state.outdir.glob("*"):
            os.remove(f)

        self.process_ligand()

        if self.form.comboBox_7.currentText() != "None":
            is_cofactor = True
            self.process_ligand(1)
        else:
            is_cofactor = False

        if self.form.comboBox_2.currentText() == "None":
            self.process_receptor()
            receptor = self.form.comboBox_1.currentText()
            gpf_out = receptor + "_preview"
            receptor = receptor + "_preview.pdbqt"
        else:
            self.flexible_receptor()
            receptor = self.form.comboBox_1.currentText()
            gpf_out = receptor + "_receptor_rigid"
            receptor_rigid = receptor + "_receptor_rigid.pdbqt"
            receptor_flex = receptor + "_receptor_flex.pdbqt"

        ligand = self.form.comboBox_3.currentText()
        ligand = ligand + "_preview.pdbqt"

        if is_cofactor:
            cofactor = self.form.comboBox_7.currentText()
            cofactor = cofactor + "_preview.pdbqt"
        else:
            cofactor = ""

        force_field = self.form.comboBox_5.currentText()

        if force_field == "ad4":
            curr_v = cmd.get_view()
            self.prepare_FF()
            out_name = self.form.comboBox_3.currentText() + "_ad4_out.pdbqt"
            cmd.delete(self.form.comboBox_3.currentText() + "_ad4_out")
            flags = self.form.lineEdit_8.text().strip()

            if self.form.comboBox_2.currentText() == "None":
                task = (
                    f"./../../vina_1.2.3_linux_x86_64 --ligand {ligand} {cofactor} "
                    f"--maps {gpf_out} --scoring ad4 --cpu {self.cpu_cnt} "
                    f"{flags} --out {out_name}"
                )
            else:
                task = (
                    f"./../../vina_1.2.3_linux_x86_64 --ligand {ligand} {cofactor}--flex {receptor_flex}"
                    f" --maps {gpf_out} --scoring ad4 --cpu {self.cpu_cnt} "
                    f"{flags} --out {out_name}"
                )

            self.sync_to_remote_and_task(task, sc=1)

            if self.is_water:
                dry_script = self.helper_scripts + "dry.py"
                task = f"{self.python3} {dry_script} -r {gpf_out + '.pdbqt'} -m {gpf_out}.W.map -i {out_name}"
                self.sync_to_remote_and_task(task)

            cmd.load(self.state.outdir / out_name)
            out_sticks = ".".join(out_name.split(".")[:-1])
            cmd.show("licorice", out_sticks)
            cmd.set_view(curr_v)

        else:
            curr_v = cmd.get_view()
            out_name = self.form.comboBox_3.currentText() + "_vina_out.pdbqt"
            cmd.delete(self.form.comboBox_3.currentText() + "_vina_out")

            flags = self.form.lineEdit_8.text().strip()
            self.create_config_file()
            if self.form.comboBox_2.currentText() == "None":
                task = (
                    f"./../../vina_1.2.3_linux_x86_64 --receptor {receptor} --ligand {ligand} {cofactor}"
                    f" --config config.txt --cpu {self.cpu_cnt}"
                    f" {flags} --out {out_name}"
                )
            else:
                task = (
                    f"./../../vina_1.2.3_linux_x86_64 --receptor {receptor_rigid} --ligand {ligand} {cofactor}"
                    f" --flex {receptor_flex} --config config.txt --cpu {self.cpu_cnt}"
                    f" {flags} --out {out_name}"
                )

            self.sync_to_remote_and_task(task, sc=1)
            self.convert_output_to_mol2(out_name)

            cmd.load(self.state.outdir / out_name)
            out_sticks = ".".join(out_name.split(".")[:-1])
            cmd.show("licorice", out_sticks)
            cmd.set_view(curr_v)

        self.save_()

    def convert_output_to_mol2(self, out_name):
        task = (
            "/bin/bash /home/ubuntu/ADFRsuite_x86_64Linux_1.0/bin/obabel -i"
            f" pdbqt {out_name} -o mol2 -O converted_poses.mol2"
        )
        self.sync_to_remote_and_task(task, [out_name], sc=1)
        cmd.delete("converted_poses")
        cmd.load(self.state.outdir / "converted_poses.mol2")
        cmd.show("licorice", "converted_poses")
        cmd.disable("converted_poses")

    def sync_to_remote_and_task(self, task, items_=None, sc=0):
        self.log(task)

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.load_system_host_keys()
        ssh.connect(
            hostname=self.hostname,
            port=self.connect_on_port,
            username=self.username,
            password=self.password,
        )

        temp_folder = time.ctime().replace(" ", "").replace(":", "") + str(
            random.randint(0, 10e8)
        )

        remote_file_path = "vina/tempfolders/" + temp_folder

        (stdin, stdout, stderr) = ssh.exec_command(command=f"mkdir {remote_file_path}")
        (stdin, stdout, stderr) = ssh.exec_command(f"chmod 0777 {remote_file_path}")
        ftp_client = ssh.open_sftp()

        if items_ is None:
            items_ = os.listdir(self.state.outdir)

        for f in items_:
            ff = os.path.join(self.state.outdir, f)
            ff_remote = remote_file_path + "/" + f
            ftp_client.put(ff, ff_remote)

        ftp_client.close()

        (stdin, stdout, stderr) = ssh.exec_command(
            command=f"cd {remote_file_path} && {task}"
        )

        for f in stdout.readlines():
            self.log(f)
            if self.verbose or sc:
                print(f)

        for f in stderr.readlines():
            self.log(f)
            print(f)

        (stdin, stdout, stderr) = ssh.exec_command(command=f"ls {remote_file_path}")

        retrieve = [x.replace("\n", "") for x in stdout.readlines()]
        retrieve = [x for x in retrieve if x not in os.listdir(self.state.outdir)]

        ftp_client = ssh.open_sftp()
        for r in retrieve:
            ff = os.path.join(self.state.outdir, r)
            ftp_client.get(remote_file_path + "/" + r, ff)
        ftp_client.close()

        (stdin, stdout, stderr) = ssh.exec_command(command=f"rm -rf {remote_file_path}")

    def hydrated_force_field(self):
        current_text = self.form.lineEdit_3.text().strip()

        if "-w" in current_text or "--hydrate" in current_text:
            self.form.comboBox_5.clear()
            self.form.comboBox_5.addItems(["ad4"])

        else:
            self.form.comboBox_5.clear()
            self.form.comboBox_5.addItems(["vina", "ad4"])

    def update_combo(self):
        selections = ["None"] + cmd.get_names("all")
        selections = [x for x in selections if not x.startswith("_")]

        try:
            selections.remove("Box")
        except ValueError:
            pass

        AllItems = [
            self.form.comboBox_1.itemText(i)
            for i in range(self.form.comboBox_1.count())
        ]

        if selections != AllItems:
            self.update_selections(self.form.comboBox_1, selections)
            self.update_selections(self.form.comboBox_2, selections)
            self.update_selections(self.form.comboBox_3, selections)
            self.update_selections(self.form.comboBox_4, selections)
            self.update_selections(self.form.comboBox_7, selections)
            # self.selections = selections

    def update_selections(self, combobox, selections):
        AllItems = [combobox.itemText(i) for i in range(combobox.count())]
        curr = combobox.currentText()
        idx = [t for t, x in enumerate(selections) if x == curr]
        idx = (idx[0]) if len(idx) == 1 else 0
        combobox.clear()
        combobox.addItems(selections)
        combobox.setCurrentIndex(idx)

    def process_receptor(self):
        self.state.check_outdir(self)

        curr_v = cmd.get_view()

        selected = self.form.comboBox_1.currentText()
        outfile = selected + "_preview.pdb"
        outfile_pdbqt = selected + "_preview.pdbqt"

        if os.path.exists(self.state.outdir / outfile):
            os.remove(self.state.outdir / outfile)
        if os.path.exists(self.state.outdir / outfile_pdbqt):
            os.remove(self.state.outdir / outfile_pdbqt)

        cmd.delete("rigid_receptor_preview")

        # remove alternate conformations
        try:
            cmd.remove("not alt ''+A")
            cmd.alter("all", "alt=''")
        except:
            print("EXCEPT alt")

        flags = self.form.lineEdit_1.text().strip()

        if selected != "None":
            cmd.save(self.state.outdir / outfile, selection=selected, state=-1)

            file_pdbqt = selected + "_preview.pdbqt"
            prep_receptor_script = self.ADFRsuite + "prepare_receptor"
            task = f"{prep_receptor_script} -r {outfile} -o {file_pdbqt} {flags} "
            self.sync_to_remote_and_task(task, [outfile])
            cmd.load(self.state.outdir / outfile_pdbqt, "rigid_receptor_preview")

            cmd.set_view(curr_v)
        else:
            print("Nothing selected")

    def flexible_receptor(self):
        self.state.check_outdir(self)

        receptor_selection = self.form.comboBox_1.currentText()
        flexible_selection = self.form.comboBox_2.currentText()

        if receptor_selection == "None":
            print("select a receptor first")
            return
        elif flexible_selection == "None":
            print("select something first")
            return

        self.process_receptor()

        receptor_pdbqt = receptor_selection + "_preview.pdbqt"

        flex_space = {"flex": [], "chain_id": []}

        cmd.iterate(flexible_selection, "flex.append(resn + resi)", space=flex_space)
        cmd.iterate(flexible_selection, "chain_id.append(chain)", space=flex_space)
        chain = flex_space["chain_id"][0]

        all_the_flex = list(set(flex_space["flex"]))

        if len(all_the_flex) > 20:
            all_the_flex = all_the_flex[:20]
            print(f"too many flex residues, max nbr: 10, reducing to: {all_the_flex}")

        flex_space = "_".join(all_the_flex)

        flex_space = receptor_selection + "_preview:" + chain + ":" + flex_space

        flexreceptor_script = self.helper_scripts + "prepare_flexreceptor.py"

        flags = self.form.lineEdit_2.text().strip()

        rigid_output = receptor_selection + "_receptor_rigid.pdbqt"
        flexible_output = receptor_selection + "_receptor_flex.pdbqt"

        if os.path.exists(self.state.outdir / rigid_output):
            os.remove(self.state.outdir / rigid_output)
        if os.path.exists(self.state.outdir / flexible_output):
            os.remove(self.state.outdir / flexible_output)

        task = f"{self.pythonsh} {flexreceptor_script} -r {receptor_pdbqt} -s {flex_space} {flags} \
             -g {rigid_output} -x {flexible_output}"

        self.sync_to_remote_and_task(task, [receptor_pdbqt])

    def process_ligand(self, lig_n=0):
        self.state.check_outdir(self)

        curr_v = cmd.get_view()

        if lig_n == 0:
            selected = self.form.comboBox_3.currentText()
            preview_name = "ligand_preview"
        else:
            selected = self.form.comboBox_7.currentText()
            preview_name = "cofactor_preview"

        outfile = selected + "_preview.pdbqt"
        cmd.delete(preview_name)

        if os.path.exists(self.state.outdir / outfile):
            os.remove(self.state.outdir / outfile)

        if selected != "None":
            file_in = selected + "_preview.sdf"
            file_out = selected + "_preview.pdbqt"

            try:
                cmd.save(self.state.outdir / file_in, selection=selected, state=-1)
            except:
                # trying to fix corrupted molecule objects in pymol. this workaround seems to work
                cmd.copy(selected, selected)
                cmd.save(self.state.outdir / file_in, selection=selected, state=-1)

            ligand_script = (
                "/home/ubuntu/miniconda3/condabin/conda run mk_prepare_ligand.py"
            )

            if lig_n == 0:
                flags = self.form.lineEdit_3.text().strip()
            else:
                flags = self.form.lineEdit_10.text().strip()

            task = f"{ligand_script} -i {file_in} -o {file_out} {flags}"

            self.sync_to_remote_and_task(task, [file_in])

            cmd.load(self.state.outdir / outfile, preview_name)
            cmd.set_view(curr_v)

        else:
            print("Nothing selected")

    def create_config_file(self):
        state = cmd.count_states("grid_center")
        coords = cmd.get_coordset("grid_center", state)
        pX, pY, pZ = [x.round(1) for x in coords.astype("float32")[0]]

        if os.path.exists(self.state.outdir / "config.txt"):
            os.remove(self.state.outdir / "config.txt")

        with open(self.state.outdir / "config.txt", "a+") as f:
            f.write(f"center_x = {pX}\n")
            f.write(f"center_y = {pY}\n")
            f.write(f"center_z = {pZ}\n")
            f.write(f"size_x = {self.form.lineEdit_4.text().strip()}\n")
            f.write(f"size_y = {self.form.lineEdit_5.text().strip()}\n")
            f.write(f"size_z = {self.form.lineEdit_6.text().strip()}\n")
            f.close()

    def prepare_FF(self):
        # self.process_ligand()

        if self.form.comboBox_2.currentText() == "None":
            # self.process_receptor()
            receptor = self.form.comboBox_1.currentText()
            gpf_out = receptor + "_preview"
            receptor = receptor + "_preview.pdbqt"

        else:
            # self.flexible_receptor()
            receptor = self.form.comboBox_1.currentText()
            gpf_out = receptor + "_receptor_rigid"
            receptor = receptor + "_receptor_rigid.pdbqt"

        ligand = self.form.comboBox_3.currentText()
        ligand = ligand + "_preview.pdbqt"

        ligand_atoms, self.is_water = self.return_atom_types(self.state.outdir / ligand)

        if self.form.comboBox_7.currentText() != "None":
            cofactor = self.form.comboBox_7.currentText()
            cofactor = cofactor + "_preview.pdbqt"
            cofactor = self.state.outdir / cofactor
        else:
            cofactor = None

        ligand_atoms, self.is_water = self.return_atom_types(
            self.state.outdir / ligand, cofactor
        )

        FF_script = self.helper_scripts + "prepare_gpf.py"

        # box coords, size, points, and intervals
        state = cmd.count_states("grid_center")
        coords = cmd.get_coordset("grid_center", state)
        pX, pY, pZ = coords.astype("float32")[0].round(2)
        gridcenter = str(pX) + "," + str(pY) + "," + str(pZ)

        spacing = 0.375

        npts_x = math.ceil(int(self.form.lineEdit_4.text().strip()) / spacing)
        npts_y = math.ceil(int(self.form.lineEdit_5.text().strip()) / spacing)
        npts_z = math.ceil(int(self.form.lineEdit_6.text().strip()) / spacing)

        npts = str(npts_x) + "," + str(npts_y) + "," + str(npts_z)

        task = f"{self.pythonsh} {FF_script} -l {ligand} -r {receptor} -v -p gridcenter={gridcenter} -p npts={npts} -p ligand_types={ligand_atoms}"
        self.sync_to_remote_and_task(task, [ligand, receptor])

        # ad4_loc = "/home/david/Desktop/vina/ADFRsuite_x86_64Linux_1.0/bin/autogrid4"
        task = f"{self.ad4_loc} -p {gpf_out + '.gpf'} -l {gpf_out + '.glg'}"
        self.sync_to_remote_and_task(task)

        if self.is_water:
            water_script = self.helper_scripts + "mapwater.py"
            out_rec = ".".join(receptor.split(".")[:-1])
            task = f"{self.python3} {water_script} -r {receptor} -s {out_rec}.W.map"
            self.sync_to_remote_and_task(task)

    def return_atom_types(self, filename, file_b=None):
        with open(filename, "r") as f:
            file = f.read()
        atom_types = list(
            set(
                [
                    f.strip().split(" ")[-1] if f.strip().startswith("ATOM") else ""
                    for f in file.split("\n")
                ]
            )
        )

        if file_b:
            with open(file_b, "r") as f:
                file = f.read()
                atom_types += list(
                    set(
                        [
                            f.strip().split(" ")[-1]
                            if f.strip().startswith("ATOM")
                            else ""
                            for f in file.split("\n")
                        ]
                    )
                )
                atom_types = list(set(atom_types))

        # make sure OA and HD are in atom types, as they define water
        if "W" not in atom_types:
            atom_types = [x for x in atom_types if x not in [""]]
            return ",".join(atom_types), 0

        else:
            atom_types = [x for x in atom_types if x not in ["", "W"]]
            atom_types += ["HD", "OA"]
            atom_types = list(set(atom_types))
            return ",".join(atom_types), 1

    def box_me(self):
        selection = self.form.comboBox_4.currentText().strip()
        if selection == "None":
            # print("cant center nothing")
            return

        if selection != "grid_center":
            cmd.delete("grid_center")
            pX, pY, pZ = cmd.centerofmass(selection)
            cmd.pseudoatom("grid_center", pos=[pX, pY, pZ])
        else:
            state = cmd.count_states("grid_center")
            coords = cmd.get_coordset(selection, state)
            pX, pY, pZ = coords.astype("float32")[0]

        self.pX = pX
        self.pY = pY
        self.pZ = pZ

        cmd.show("spheres", "grid_center")

        self.update_box()

    def update_box(self):

        if not hasattr(self, "pX"):
            return

        curr_v = cmd.get_view()
        cmd.delete("Box")

        x_temp = self.form.lineEdit_4.text().strip()
        y_temp = self.form.lineEdit_5.text().strip()
        z_temp = self.form.lineEdit_6.text().strip()

        x = self.x if x_temp == "" else x_temp
        self.x = x
        y = self.y if y_temp == "" else y_temp
        self.y = y
        z = self.z if z_temp == "" else z_temp
        self.z = z

        minX = self.pX + float(x)
        minY = self.pY + float(y)
        minZ = self.pZ + float(z)
        maxX = self.pX - float(x)
        maxY = self.pY - float(y)
        maxZ = self.pZ - float(z)

        boundingBox = [
            BEGIN,
            LINES,
            VERTEX,
            minX,
            minY,
            minZ,
            VERTEX,
            minX,
            minY,
            maxZ,
            VERTEX,
            minX,
            maxY,
            minZ,
            VERTEX,
            minX,
            maxY,
            maxZ,
            VERTEX,
            maxX,
            minY,
            minZ,
            VERTEX,
            maxX,
            minY,
            maxZ,
            VERTEX,
            maxX,
            maxY,
            minZ,
            VERTEX,
            maxX,
            maxY,
            maxZ,
            VERTEX,
            minX,
            minY,
            minZ,
            VERTEX,
            maxX,
            minY,
            minZ,
            VERTEX,
            minX,
            maxY,
            minZ,
            VERTEX,
            maxX,
            maxY,
            minZ,
            VERTEX,
            minX,
            maxY,
            maxZ,
            VERTEX,
            maxX,
            maxY,
            maxZ,
            VERTEX,
            minX,
            minY,
            maxZ,
            VERTEX,
            maxX,
            minY,
            maxZ,
            VERTEX,
            minX,
            minY,
            minZ,
            VERTEX,
            minX,
            maxY,
            minZ,
            VERTEX,
            maxX,
            minY,
            minZ,
            VERTEX,
            maxX,
            maxY,
            minZ,
            VERTEX,
            minX,
            minY,
            maxZ,
            VERTEX,
            minX,
            maxY,
            maxZ,
            VERTEX,
            maxX,
            minY,
            maxZ,
            VERTEX,
            maxX,
            maxY,
            maxZ,
            END,
        ]

        boxName = "Box"
        cmd.load_cgo(boundingBox, boxName)
        cmd.set_view(curr_v)

    def help_receptor(self, info_file):
        global window_
        window_ = QtWidgets.QDialog()
        uifile_ = os.path.join(os.path.dirname(__file__), "ui_files")
        uifile_ = os.path.join(uifile_, info_file)
        form_ = loadUi(uifile_, window_)
        window_.show()

    def mousePressEvent2(self, event):
        return

    # https://stackoverflow.com/questions/63667232/how-to-move-window-when-dragging-frame-pyqt5
    def mousePressEvent(self, event):
        self.dragPos = event.globalPos()

    def mouseMoveEvent(self, event):
        if event.buttons() == QtCore.Qt.LeftButton:
            self.move(self.pos() + event.globalPos() - self.dragPos)
            self.dragPos = event.globalPos()
            event.accept()
