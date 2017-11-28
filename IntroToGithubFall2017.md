![HMS](./img/logo.png)

# HMS Research Computing: Introduction to Git and Github

## Welcome

**Term: Fall 2017**

**Instructor: Mike McClellan <mike_mcclellan@hms.harvard.edu>**

**Date: 2017-11-29 3pm-5pm**

Welcome to Introduction to Git and Github. This is our first time offering this course and we hope to improve the course over time. At the end of the session there will be a link to a course survey; any feedback you can offer would be greatly appreciated.

## Registering for Github

First things first, we should probably register for a Github account if you don't already have one. Registering for Github is free and easy. All we need is a username, email, and password. Github sends a **Verification Email** that may take a couple of minutes, so we should probably get that process started.

Please navigate to [Join GitHub](https://github.com/join). The page should look like the figure below:

![Register for Github](./img/hmsrcght-register.png)

There are also paid plans for Github described at [Plans for all workflows](https://github.com/pricing). A paid Developer account is $7 per month and the main feature is unlimited private repositories. Team accounts start at $25 for 5 users and in addition to unlimited private repositories, you also get team and user permissions, which makes access control easy to setup. That is the way we manage access to our repositories in [Research Computing](https://github.com/hmsrc).

## What is Git

So I started preparing for this class the way I normally research a new topic:

![What is Git](./img/hmsrcght-whatis.png)

Hmm. Ok. Well the Git we will be talking about today is actually an exciting piece of software; and its the most widely used modern version control system in the world.   

It was originally developed by Linus Torvalds. You may have heard of his other popular project, the [Linux kernel](https://github.com/torvalds/linux) (which incindentally can be found on Github).

Git is actually a Distributed Version Control System. The key word being *distributed*. If you've ever worked on older version control systems, like CVS or Sourcesafe, you might remember that these systems were centralized and had one master place for the full version history of a repository. With Git, every collaborator has a repository that can contain the full history of all changes.

The details of how Git came to be are interesting. The Linux project was actually using a proprietary version control system and a licensing change occurred. Details are available in a number of places including the [Wikipedia](https://en.wikipedia.org/wiki/Git) article. There are also many places to learn more about Git including it's [main site](https://git-scm.com/).

Interestingly, if you go to the man page for Git (by typing `man git`), you will see that it refers to itself as: **git - the stupid content tracker**

The man page also informs us that *Git is a fast, scalable, distributed revision control system with an unusually rich command set*.

And that is true. But it is also true that learning just a dozen simple commands will provide you with 90%+ of everything you'll ever want to do in Git.

We'll cover each of those commands and provide some additional suggestions about useful flags and options.

## But why even Version Control?

To avoid this:

![Final.doc](./img/hmsrcght-final.gif)

But seriously, version control is especially important for Research Computing. It helps you produce better, more reproducible research and makes it easier for others to collaborate with you. From *Best Practices for Scientific Computing*:

> In practice, everything that has been created manually
should be put in version control (5.2), including programs,
original field observations, and the source files for papers.
Automated output and intermediate files can be regenerated
at need.

I highly recommend the more practical version of this paper, [Good Enough Practices for Scientific Computing](http://swcarpentry.github.io/good-enough-practices-in-scientific-computing/). It's awesome. It goes in to more detail than just version control, with information on data management, tidy data, etc; but again it's a very good resource for anyone doing research computing.

## What is Github

[Github](https://github.com/) is a site that makes collaborating on software projects easier than ever. It hosts Git repositories and provides a nice web-based graphical user interface. It goes beyond that by adding useful project management features like access control, and issue tracking. It also provides tools for collaborators to make pull requests, initiate code reviews and merge changes.

We will cover some of the most useful features of Github collaboration.

## Download GitHub Desktop

It is important to note that you can do all your work directly from the `git` command line interface (CLI). Taking this further, you could even do all your work directly from Orchestra.

But realistically that is not how most newer users choose to work. A common practice that using a system like Git allows is to development locally, using some mix of the command line and graphical tools and editors and then pushing changes up to Guthub.

With that in mind, Github has developed  [GitHub Desktop](https://desktop.github.com/), a GUI version available for both **Mac** and **Windows**. I recommend trying it and seeing if it works for you.

![Github Desktop](./img/hmsrcght-desktop.png)

I will be doing demonstrations today with Github Desktop for Mac. Additionally I will use a free, open source text editor named [Atom](https://atom.io/), which is also developed by Github and integrates really well with Git. That being said, I very often work directly with the CLI  and I will attempt to show both sides.

On Linux, the `git` command is usually installed by default, as it is on Orchestra. To check if you already have Git by typing `git --version`. On Redhat/CentOS, you can install via `sudo yum install git` or on Ubuntu/Debian via `sudo apt-get install git`.

There are many additional clients and helper tools available, such as fancy prompts and fancy diffs, and you can find more info via curated lists like [awesome-git](https://github.com/dictcp/awesome-git).

## Git Configuration

Once you have registered and verified your new account, you are ready to configure git. Like everything else in Git, there are many configuration options that can be set, however we need only a few to get going.

The main thing to be aware of is that there are three places that Git checks for configurations:

- `/etc/gitconfig` - This is for system level configuration. Orchestra has many users, and there is no default system level Git configuration.

- `~/.gitconfig` - This is the global configuration and applies to all of your repositories. This is what most people want.

- `<repo>/.git/config` - This is local, per repository configuration.

The scoping works in such a way that local settings will override global settings, which override any system settings.

So to configure your name and email in the terminal for all of your Git repositories, just type:

```bash
$ git config --global user.name mmcclellan
$ git config --global user.email mcclellan.m@gmail.com
$ cat ~/.gitconfig
[user]
	name = mmcclellan
	email = mcclellan.m@gmail.com
```

There are other useful settings we may wish to configure. One important settings, especially if you are running on Windows, is to automate your line endings. This eliminates frustating differences between CR, CRLF and LF line endings. We also want to turn on color output and maybe add our favorite editor:

```bash
$ git config --global core.autocrlf true
$ git config --global color.ui auto
$ git config --global core.editor nano
$ cat ~/.gitconfig
[user]
	name = mmcclellan
	email = mcclellan.m@gmail.com
[core]
	editor = nano
	autocrlf = true
[color]
	ui = auto
```

These will get us started, but remember that there are  many additional settings you can add as you go. One popular setting is to [Cache your Github password](https://help.github.com/articles/caching-your-github-password-in-git/#platform-all). This recent article: https://blog.scottnonnenberg.com/better-git-configuration/ has some good suggestions for further configuration.

It is also very easy to do this in Github Desktop. On the Mac version you want to go to `Preferences`,  add your Github accound and set your git config:

![Github Desktop Preferences](./img/hmsrcght-config.png)

## Creating a New Repository

The first and simplist way to create a new repository is with `git init`.

Normally, `git init` is a one time operation. You may want to keep all your repositories organized under one directory. To do this start by creating a directory; maybe `projects` or `src`:

```bash
$ mkdir -p projects
$ cd projects/
$ git init project2/
Initialized empty Git repository in /Users/mmcclellan/src/project2/.git/
$ cd project2
```

Git will create a directory for your new project and inside it add a hidden .git directory:

```bash
$ ls -la
total 38
drwxrwsr-x 3 mm713 mm713  22 Apr  9 13:14 .
drwxrwsr-x 3 mm713 mm713  26 Apr  9 13:14 ..
drwxrwsr-x 6 mm713 mm713 167 Apr  9 13:14 .git
```

Feel free to inspect the contents of your .git folder with a command like `tree .git`. This folder holds ALL of the repository information. If you even wanted to "ungit" a repository, you could just `rm -rf .git/` (but of course be very careful doing so).

If you already have a directory with stuff in it, you can cd into that directory and just type: `git init` omitting a name and git will initialize a new repository right there in the current directory.

As you probably expect, creating a repository in Github Desktop is also easy; just look for the **+** sign in the upper left-hand corner of the app.

Another cool way to create a new repository is right from your github.com account. I suggest you give this a try because it adds a couple nice features. Look for the **+** again, this time on the right-hand side. The first thing it prompts for is Repository Name. Your Github username creates a **namespace** so the repository name only has to be unique within your account. While description is also, I think its nice to add one. Then there is a decision to make the repository Public or Private. Private repositories are only available to Paid accounts. I find the cost well worth the investment. Then there is a checkbox for initializing a README file. This is a conventionl and a very worthwhile one in my opinion. Traditionally, it uses all CAPS so that on Unix file systems it would normmally appear first in the file listing. One of Github's founders wrote an interesting piece on [Readme Driven Development](http://tom.preston-werner.com/2010/08/23/readme-driven-development.html). You can also add a `/.gitignore` file right away (we will cover /gitignore later on). And lastly you can choose a LICENSE. Let's not get into licensing issues right now, except that to say that all public repositories should have a license, whichever one that may be.

![github.com create new repository](./img/hmsrcght-createnew.png)

## `git status` is your friend

This is just a quick aside. Seriously, just check `git status` often. It's really the fastest way to learn what state Git "thinks" your repository is in. Let's check the new repository's status:

```bash
$ git status
On branch master

Initial commit

nothing to commit (create/copy files and use "git add" to track)
```

## Adding a new file to a repository

Let's create a first file in the terminal:

```bash
$ touch list.md
$ git status
On branch master

Initial commit

Untracked files:
  (use "git add <file>..." to include in what will be committed)

	list.md
```

Git recognizes that something has changed in the directory, but as of yet, the file is untracked. Let's make a couple of edits to the file. I will use the nano editor.

```bash
$ nano list.md
$ cat list.md
# List 1

This is a first list for project2

- item 1
- item 2
- item 4
```

Remember, this file is untracked, and we can make all the edits we want without it being tracked. When we are ready to start tracking it, we do a git add:

```bash
$ git add list.md
$ git status
On branch master

Initial commit

Changes to be committed:
  (use "git rm --cached <file>..." to unstage)

	new file:   list.md
```

This is one of the very most important things about Git; that staging a file (via `git add`) is not the same as committing it. Staging a file makes it eligible to be committed. You cannot commit a file without adding (staging) it.

The Github Desktop GUIs actually do allow to add and commit in one step, but behind the scenes the add and commit are separate.

## Committing our first file

The file is now staged for committing. To make our first commit from the terminal, we just need to use the commit flag and include a commit message. The easiest way to include a message is with the `-m` flag:

```bash
$ git commit -m "My first commit to project2"
[master (root-commit) 166b9f9] My first commit to project2
 1 file changed, 7 insertions(+)
 create mode 100644 list.md
$ git status
On branch master
nothing to commit, working tree clean
```

This shows our commit and that we now do not have any uncommitted changes.

## Git Log

At any point we can review the history of the repository by using `git log` from the terminal:

```bash
$ git log
commit 166b9f9098911e721f11ccb9dde280fecf067ead
Author: Michael McClellan <mcclellan.m@gmail.com>
Date:   Sun Apr 9 19:45:44 2017 -0400

    My first commit to project2
```

The log shows a commit hash that uniquely identifies the commit. It also shows the Date and Author as well as our commit message.

Later, we will see that the log of commits is also formatted and displayed nicely at github.com. You can click on a commit and see the full diff and all the metadata about the commit.

### File Modifications

Well it looks like we have our first bug. I omitted item 3 from the list. After editing, I can see what the difference between my working tree and the last commit by using the `git diff` command:

```bash
$ nano list.md
$ git diff list.md
diff --git a/list.md b/list.md
index 715c9de..9a4597a 100644
--- a/list.md
+++ b/list.md
@@ -4,4 +4,5 @@ This is a first list for project2

 - item 1
 - item 2
+- item 3
 - item 4
```

It shows that I added a line to the file. You can also use the `--stat` command to suppress the patch information. Let's also check status again:

```bash
$ git status
On branch master
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   list.md

no changes added to commit (use "git add" and/or "git commit -a")
```

Git recognizes that something has been changed. This time its not a new file but a modification of an existing file. The new change has not been staged though, so to stage it we again use `git add`. It may seem weird to use add again, but this is how git differentiates untracked changes from those that are staged and ready to be committed.

```bash
$ git add list.md
$ git commit -m "Fixed missing item"
[master 7bc5ee0] Fixed missing item
 1 file changed, 1 insertion(+)
```

In a real project you may have many files in different states with some ready to be committed and others not ready. Let's make a few more edits for an example scenario:

```bash
$ nano list2.md
$ nano list.md
$ git diff list.md
diff --git a/list.md b/list.md
index 9a4597a..1bbe450 100644
--- a/list.md
+++ b/list.md
@@ -6,3 +6,4 @@ This is a first list for project2
 - item 2
 - item 3
 - item 4
+- item Im still thinking about
$ git status
On branch master
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   list.md

Untracked files:
  (use "git add <file>..." to include in what will be committed)

	list2.md

no changes added to commit (use "git add" and/or "git commit -a")
```

So see Git sees both changes, but maybe I'm only ready to commit list2.md. In that case I just add it and commit it without re-staging list.md:

```bash
$ git add list2.md
$ git commit -m "Adding second list"
[master c0584c2] Adding second list
 1 file changed, 8 insertions(+)
 create mode 100644 list2.md
$ git status
On branch master
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   list.md

no changes added to commit (use "git add" and/or "git commit -a")
```

Let's say I decided the addition to list.md was a mistake. I can discard the difference like so:

```bash
$ git checkout -- list.md
$ git diff list.md
$ git status
On branch master
nothing to commit, working tree clean
```

Let's check the log again:

```bash
$ git log
commit c0584c2e22b44abc0e230b6d032ad55d9b606cf7
Author: Michael McClellan <mcclellan.m@gmail.com>
Date:   Sun Apr 9 20:02:12 2017 -0400

    Adding second list

commit 7bc5ee0cc0981e2b20e0ac396a84e5a6bf73c092
Author: Michael McClellan <mcclellan.m@gmail.com>
Date:   Sun Apr 9 19:52:46 2017 -0400

    Fixed missing item

commit 166b9f9098911e721f11ccb9dde280fecf067ead
Author: Michael McClellan <mcclellan.m@gmail.com>
Date:   Sun Apr 9 19:45:44 2017 -0400

    My first commit to project2
```

I now have three commits and they are displayed in chronological order.

## Working with Github Desktop

So we could happily continue on - working in the terminal forever, but there are some features of the GUI client that are really nice to have. So we have decided to make our repository managable in Github Desktop.

Open Github Desktop and choose File > Add Local Repository. Navigate to the repository you created and select it. In the left-hand pane, is a list of repos you have access to on Github and below that is Other, where your local repo should show up (because it is not up on github.com yet).

We should also be able to see the history and each commits change set as well. Navigate back to the terminal and add a couple of items to list2.md like so:

```bash
$ nano list2.md
$ git diff list2.md
diff --git a/list2.md b/list2.md
index bde4d6a..dee3a05 100644
--- a/list2.md
+++ b/list2.md
@@ -6,3 +6,5 @@ Here are some items we need:
 - butter
 - sugar
 - coffee
+- flour
+- salt
```

Then toggle back to Github Deskop and you should see one uncommitted change. We can commit in the GUI in one step by adding a commit message and clicking commit to master. It should now also show up in the history pane.

## Publish to Github

On the right side of the GUI, you should see a **Publish** button. When we click it a dialog box should pop-up with Repository Name and description boxes, a private repository checkbox, and a drop-down Account menu (if you have write access to more than one Github account), Fill out the details and click Push to Github.

Back in the left-hand pane we should see the repo has now moved from the Other section to the Github section. Right click on the repository and select **View on Github**.

![view on github.com](./img/hmsrcght-openongithub.png)

A browser should open with your new Github repository. Take a minute to review your commits and look around the repo.

## Github Flow

As teams starting making Git a core part of their development process, several Git workflows became popular. One in particular, [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/), gained a lot of attention, despite its complexity:

![view on github.com](./img/hmsrcght-gitflow.png)

Well. OK.

So, [Github Flow](https://guides.github.com/introduction/flow/) emerged as a simplified alternative to Git Flow.

And simple it is. Now that our repo is up on Github, let's utilize Github Flow for our project.

## Github Issues

You've probably use some issue tracker before like Jira, or, ahem **Stat**. Well Github Issues is like that, only simpler.

Let's create a new issue for some enhancement. In your repo at github.com, click Issues and then click **New Issue**.

![New Issue](./img/hmsrcght-newissue.png)

Fill out the title and description. If you want, assign it to yourself and choose a label like **enhancement**. Then click **Submit New Issue**.

## Creating a New Branch

Now that we have an Issue assigned to us, rather than developing right on the `master` branch, let's create a new branch.

We can of course create a branch from the terminal with `git branch -d <branch>` but let's try it inside of Github Desktop. In the top pane next to the branch drop-down (which should only have master at this point) , there is a branch button that allows you to create a new branch. Let's create a new branch called `add_todo_list`. Now toggle back to the terminal and check status:

```bash
$ git status
On branch add_todo_list
nothing to commit, working tree clean
```

We have been moved to our new branch. You can always see your branches like so:

```bash
$ git branch
* add_todo_list
  master
```

And toggle between branches using checkout:

```bash
$ git checkout master
Switched to branch 'master'
Your branch is up-to-date with 'origin/master'.
$ git branch
  add_todo_list
* master
$ git checkout add_todo_list
Switched to branch 'add_todo_list'
$ git branch
* add_todo_list
  master
```

Open `todo.md` and make some edits and save. If we switch back to Github Desktop, we can see our uncommitted changes. Let's do the commit from Github Desktop. Now we have a commit on our branch.

At this point, the branch is only local, it has not been published to Github yet. We can do so by clicking **Publish** in Github Desktop.

Wait a few seconds and then View the Repository on Github. It should show our new branch as recently published and add a button for **Compare and Pull Request**.

![Create PR](./img/hmsrcght-createpr.png)

Clicking that button should open a Pull Request form, where you should type `Fixes #`, which should pull up a list of issues. Select the issue we added and then click the **Create Pull Request** button.

Github will make sure that there are no conflicts. In a real project, we would want a collaborator to review our Pull Request. You can even force code reviews in the repository settings before a merge can be made.

In our case let's go ahead and merge. After the merge, it offers to delete the branch and if you are done with it, go ahead and do so.

![Merge PR](./img/hmsrcght-mergepr.png)

If we Navigate back to Issues, we see that the issue was closed for us. This was done by the Fixes #1 comment in our Pull Request.

![Closed Issue](./img/hmsrcght-closedissue.png)

These are the basic steps of the Github Flow process. Create an Issue, Create a Branch, Make Needful Commits, when satisfied create a Pull Request, Review the Changes, Merge the changes, and close the Issue.

While the are other features and capabilities in Github, just these simple steps go along ways towards establishing a collaborative coding environment.

## Cloning an existing repository

It is also easy to clone an existing repository from either the terminal or Github Desktop. First, navigate to https://github.com/hmsrc/hmsrcght.

This is a sample repo I set up for this exercise. Look for a green button that says **Clone or download**. To clone via Github Desktop, just click Open in Desktop. To clone from the terminal, copy the repo URL and then in the terminal, type: `git clone https://github.com/hmsrc/hmsrcght.git`.

You may have noticed that Github also supports Forking a Repository, which is like cloning to your local workstation, but makes Github.com "aware" of your server side clone.

I wanted to point out a couple of nice rendering features of Github, like support for Jupyter Notebooks.

I also wanted to add a few comments about `.gitignore`.

## `gitignore` and What Not to Share

It is really important to note that researchers dealing with data subject to legal restrictions that prohibit sharing (such as medical data) should be careful not to put sensitive data into repositories.

Also, be extra careful not to accidentally put credentials, such as passwords and private keys, into repositories.

Additionally, intermediate data files and other results that can be re-generated from raw data need not be added to version control.

Adding a .gitignore file to your repo is a good way to help ensure you don't inadvertantly share sensitive info.

## Github Organizations

So I wanted to quickly show you how a Github Organization like [HMS Resource Computing](https://github.com/hmsrc) can really improve collaboration.

You can organize your collaborators into teams and then assign whole teams granular access to your repositories. You can grant read, write, and admin priviledges to each repository.

![Teams](./img/hmsrcght-teams.png)

Many of our repositories have thousands of commits. You can track issues, as we've seen, contributions:

![Contrib](./img/hmsrcght-contrib.png)

and many other things. You can use Trello-style Kanban boards for light weight project management:

![Kanban](./img/hmsrcght-kanban.png)

I feel like I'm forgetting an important collaboration feature.

## Oh yeah, Emoji :+1:

:bowtie: | :smile: | :laughing: | :blush: | :smiley: | :relaxed: | :smirk: | :heart_eyes: | :flushed: | :relieved: | :satisfied: | :grin: | :wink: |  :tongue: | :unamused: | :sweat_smile: | :sweat:


:+1: | :thumbsdown: | :ok_hand: | :punch: | :facepunch: | :fist: | :v: | :wave: | :hand: | :open_hands: | :point_up: | :point_down: |  :raised_hands: | :pray: | :point_up_2: | :clap: | :muscle: | :metal:

## Q & A

I wanted to leave some time open for questions and answers. For questions we can't answer right away, we will capture and follow up.

## Contact information

- [Orchestra New User Guide](https://wiki.med.harvard.edu/Orchestra/NewUserGuide)
- [HMS RC website](http://rc.hms.harvard.edu)
- [rchelp@hms.harvard.edu](mailto:rchelp@hms.harvard.edu)
- Office Hours: Wednesdays 1-3p Gordon Hall 500
- [Class Survey](http://hmsrc.me/introgithub2017-survey1)
